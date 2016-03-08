package demo;

import java.util.ArrayList;
import java.util.List;

import psimj.Communicator;
import psimj.LocalCommunicator;
import psimj.Task;
import psimj.Topology;
import psimj.Topology.TopologyViolationException;

public class Demo2 {

	public static class MyTask implements Task {

		@SuppressWarnings("unchecked")
		@Override
		public void run(Communicator comm) {
			try {
				
				// BROADCAST
				// Rank 0 populates a Java object with random numbers
				ArrayList<Integer> list = null;
				if (comm.rank() == 0) {
					list = new ArrayList<Integer>();
					for (int i = 0; i < 10; i ++) {
						list.add((int)(Math.random() * 100));
					}
				}
				
				// Rank 0 broadcasts the Java object, all Ranks receive it
				list = comm.one2all_broadcast(0, list, ArrayList.class);
				
				// Format a string using those numbers
				String numbers = "";
				for (Integer i : list) {
					numbers += i + " ";
				}
				
				// Print the received numbers
				System.out.println("Rank " + comm.rank() + " got: " + numbers);
				
				// COLLECT
				// All Ranks send their integer Rank to Rank 0
				List<Integer> ranks = comm.all2one_collect(0, comm.rank(), int.class);
				
				// Rank 0 prints
				if (comm.rank() == 0) {
					numbers = "";
					for (Integer i : ranks) {
						numbers += i + " ";
					}
					System.out.println("\nRank " + comm.rank() + " collected: " + numbers);
				}
			} catch (TopologyViolationException e) {
				e.printStackTrace();
			}
			comm.finish();
		}

	}

	public static void main(String[] args) throws Exception {
		
		/*
		 * Demo 2
		 * 
		 * Using broadcast and collect.
		 * 
		 */
		
		Communicator comm = new LocalCommunicator(12, new Topology.Switch());
		
		comm.runTask(MyTask.class);
		
	}
}
