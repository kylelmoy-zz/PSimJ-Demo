package demo;

import java.net.InetAddress;
import java.util.ArrayList;
import java.util.List;

import psimj.Communicator;
import psimj.Task;
import psimj.network.NetworkCommunicator;
import psimjpool.Pool;
import psimjpool.PoolKey;

public class Demo3 {

	public static class MyTask implements Task {

		@Override
		public void run(Communicator comm) {
			try {
				// Find out this computer's name
				String myName = InetAddress.getLocalHost().getHostName();
				
				// Collect all names to Node 0
				List<String> names = comm.all2one_collect(0, myName, String.class);
				
				// Node 0 prints all names
				if (comm.rank() == 0) {
					System.out.println("Connected nodes:");
					for (String s : names) {
						System.out.println(s);
					}
				}
			} catch (Exception e) {}
			comm.finish();
		}

	}

	public static void main(String[] args) throws Exception {
		
		/*
		 * Demo 3
		 * 
		 * Using the NetworkCommunicator type Communicator, facilitated by PSimJ Pool
		 * 
		 */
		
		// Connect to a PSimJ Pool
		Pool pool = new Pool("127.0.0.1", 8191);
		
		// Authenticate using a key
		pool.useAuthentication(new PoolKey("./pseudo.key"));
		
		// Manually initiate communications
		Communicator comm = new NetworkCommunicator(new ArrayList<String>(), 0);
		
		// Use the connected Pool to facilitate constructing a network of NetworkCommunicators
		// Request 4 nodes from the Pool
		// Communicator comm = pool.requestCommunicator(4);
		
		// Submit the Task if request is successful
		if (comm != null) {
			comm.runTask(MyTask.class);
		}
		
	}
}
