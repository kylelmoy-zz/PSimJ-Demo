package demo;

import psimj.Communicator;
import psimj.LocalCommunicator;
import psimj.Task;
import psimj.Topology;
import psimj.Topology.TopologyViolationException;

public class Demo1 {

	public static class MyTask implements Task {

		@Override
		public void run(Communicator comm) {
			try {
				//Everyone sends Rank 0 a message
				comm.send(0, "Hello from rank " + comm.rank() + "!");
				
				if (comm.rank() == 0) {
					// Rank 0 receives a message from all Ranks
					for (int i = 0; i < comm.nprocs(); i++) {
						String message = comm.recv(i, String.class);
						System.out.println("Received from rank " + i + ": " + message);
					}
				}
			} catch (TopologyViolationException e) {}
			comm.finish();
		}

	}

	public static void main(String[] args) throws Exception {
		
		/*
		 * Demo 1
		 * 
		 * Using the communicator to send and receive objects.
		 * 
		 */
		
		Communicator comm = new LocalCommunicator(12, new Topology.Switch());
		
		comm.runTask(MyTask.class);
		
	}
}
