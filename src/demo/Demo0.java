package demo;

import psimj.Communicator;
import psimj.LocalCommunicator;
import psimj.Task;
import psimj.Topology;

public class Demo0 {

	public static class MyTask implements Task {

		@Override
		public void run(Communicator comm) {
			// Every Rank prints a message
			System.out.println("Hello from rank " + comm.rank() + "!");
			
			comm.finish();
		}

	}

	public static void main(String[] args) throws Exception {
		
		/*
		 * Demo 0
		 * 
		 * Use a LocalCommunicator type Communicator to simulate machines using threads
		 * 
		 */
		
		// Construct a LocalCommunicator using (12) threads, and the (Switch) topology
		Communicator comm = new LocalCommunicator(12, new Topology.Switch());
		
		// Run a Task using the Communicator
		comm.runTask(MyTask.class);
		
	}
}
