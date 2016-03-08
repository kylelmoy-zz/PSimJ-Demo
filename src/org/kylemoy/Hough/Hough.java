package org.kylemoy.Hough;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.imageio.ImageIO;

import psimj.Communicator;
import psimj.LocalCommunicator;
import psimj.Task;
import psimj.Topology;
import psimjpool.Pool;
import psimjpool.PoolKey;

public class Hough implements Task {
	
	static String path = "input512.jpg";
	
	public static void main(String[] args) throws Exception {
		Pool pool = new Pool("127.0.0.1", 8191);
		pool.useAuthentication(new PoolKey("./pseudo.key"));
		//BufferedImage img = ImageIO.read(new File(path));
		//Hough h = new Hough();

		
		Communicator comm = new LocalCommunicator(1, new Topology.Switch());
		//Communicator comm = pool.requestCommunicator(4);
		if (comm != null) {
			comm.runTask(Hough.class);
		}
		
	}

	HashMap<Integer,List<Integer[]>> circleKernel = new HashMap<Integer, List<Integer[]>>();
	
	@SuppressWarnings("unchecked")
	@Override
	public void run(Communicator comm) {
		try {
			long time = System.currentTimeMillis();
			byte[] data = null;
			if (comm.rank() == 0) {
				try {
					data = Files.readAllBytes(Paths.get(path));
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
			data = comm.one2all_broadcast(0, data, byte[].class);
			BufferedImage img = ImageIO.read(new ByteArrayInputStream(data));
			int w = img.getWidth();
			int h = img.getHeight();
			int n = Math.max(w, h);
			
			int r_max = (int)(Math.ceil(Math.sqrt(Math.pow(w, 2) + Math.pow(h, 2)))) / 2;
			int r_start = 5;
			int r_end = r_max;
			//Cheating
			//r_max = 64;
			
			if ( n > 275) {
				//Optimize distribution of work based on an empirical model
				double estimatedTotalCost = estimateCost(n, 5, r_max);
				int costPerNode = (int)Math.ceil(estimatedTotalCost / (float)comm.nprocs());
				int start = 20;
				int end;
				for (int i = 0; i < comm.nprocs() - 1; i++) {
					end = (int) Math.ceil(solveRange(n, costPerNode, start));
					if (i == comm.rank()) {
						r_start = start;
						r_end = end;
					}
					if (comm.rank() == 0) {
						System.out.println("Node " + i + " gets " + start + " to " + end + " (" + (end - start) + " total)");
					}
					start = end;
				}
				if (comm.rank() == comm.nprocs() - 1) {
					r_start = start;
					r_end = r_max;
				}
				if (comm.rank() == 0) {
					System.out.println("Node " + (comm.nprocs() - 1) + " gets " + start + " to " + r_max + " (" + (r_max - start) + " total)");
					System.out.println("Estimated worst run time: " + timeFormat(costPerNode));
					//System.out.println("Estimated worst CPU time: " + timeFormat((int)estimatedTotalCost));
				}
			} else {
				//Model is bad for smaller values of N, resort to equal distribution
				int rValuesPerNode = (int)Math.ceil((r_max / (double)comm.nprocs()));
				r_start = (comm.rank() * rValuesPerNode);
				r_end = (((comm.rank() + 1) * rValuesPerNode));
			}
			
			//Utilize cores
			//Optimally, we should use the same model, but oh well
			int cores = Runtime.getRuntime().availableProcessors();
			ExecutorService executor = Executors.newFixedThreadPool(cores);
			Runnable[] threads = new Runnable[cores];
			int r_total = r_end - r_start;
			int rValuesPerThread = (int)Math.ceil((r_total / (double)cores));
			
			ArrayList<ArrayList<Integer[]>> results = new ArrayList<ArrayList<Integer[]>>();
			for (int i = 0; i < cores; i++) {
				final int j = i;
				final int rs = r_start;
				threads[i] = new java.lang.Runnable() {
					@Override
					public void run() {
						try {
							int k = j;
							//System.out.println("Thread " + k + " gets " + (rs + (k * rValuesPerThread)) + " to " + (rs + ((k + 1) * rValuesPerThread)));
							ArrayList<Integer[]> maxima = hough(img, rs + (k * rValuesPerThread), rs + ((k + 1) * rValuesPerThread));
							//System.out.println("Thread " + k + " finds " + maxima.size() + " maxima");
							synchronized(results) {
								results.add(maxima);
							}
						} catch (Exception e) {
							e.printStackTrace();
						}
					}};
				executor.submit(threads[i]);
			}
			executor.shutdown();
			executor.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
			
			ArrayList<Integer[]> maxima = new ArrayList<Integer[]>();
			for (ArrayList<Integer[]> result : results) {
				maxima.addAll(result);
			}
			
			//Single thread
			//ArrayList<Integer[]> maxima = hough(img, r_start, r_end);
			
			
			System.out.println("Node " + comm.rank() + " complete. Time: " + timeFormat((int)(System.currentTimeMillis() - time)));
			List<ArrayList<Integer[]>> maximaList = (List<ArrayList<Integer[]>>) comm.all2one_collect(0, maxima, maxima.getClass());
			
			if (comm.rank() == 0) {
				ArrayList<Integer[]> allMaxima = new ArrayList<Integer[]>();
				for (ArrayList<Integer[]> list : maximaList) {
					for (Integer[] m : list) {
						allMaxima.add(m);
					}
				}
				System.out.println("Maxima found: " + allMaxima.size());
				File f = new File("overlay.png");
				drawOverlay(img, allMaxima, w, h, f);
				System.out.println("Overlay saved as " + f.getAbsolutePath());
				System.out.println("Total elapsed time: " + (System.currentTimeMillis() - time) + "ms");
			}
			comm.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public static String timeFormat(int millis) {
		return String.format("%02d:%02d:%02d", 
			    TimeUnit.MILLISECONDS.toHours(millis),
			    TimeUnit.MILLISECONDS.toMinutes(millis) - 
			    TimeUnit.HOURS.toMinutes(TimeUnit.MILLISECONDS.toHours(millis)),
			    TimeUnit.MILLISECONDS.toSeconds(millis) - 
			    TimeUnit.MINUTES.toSeconds(TimeUnit.MILLISECONDS.toMinutes(millis)));
	}
	/**
	 * Solves the equation in estimateCost for a specific r_end given an r_start
	 * @param n
	 * @param cost
	 * @param r_start
	 * @return
	 */
	public double solveRange(int n, int cost, int r_start) {
		// Solve the function in estimateCost() for r_end
		// See estimateCost() for legibility
		double a = 35.9225;
		double b = (-0.635500374085847) + ((0.00200790898982706) * n) + ((0.0000132198262616506) * Math.pow(n,2));
		double c = (-0.16741758251155) + (0.209443008293443) * Math.log(n);
		double y = cost;
		double z = r_start;
		
		// Don't ask
		double x = (Math.sqrt(c * ((Math.pow(a, 2) * c) + (2 * a * b * c * z) + (Math.pow(b, 2) * c * Math.pow(z, 2)) + (2 * b * y))) - (a * c)) / (b * c);
		if (x > r_start) {
			return x;
		} else {
			//Only occurs when N < ~275 ish
			return r_start + 1;
		}
	}
	/**
	 * Calculate an estimate worst case time for an n * n image and a range of Rs
	 * @param n
	 * @param r_start
	 * @param r_end
	 * @return
	 */
	public double estimateCost(int n, int r_start, int r_end) {
		// Technically this returns an exact time in milliseconds
		// But we're only concerned with relative "cost"
		
		// Empirically determined quadratic model of the relationship between N and the coefficient in the "time(r)" model
		double b0 = -0.635500374085847;
		double b1 = 0.00200790898982706;
		double b2 = 0.0000132198262616506;
		double tprCoefficient = b0 + (b1 * n) + (b2 * Math.pow(n,2));
		
		// Linear "time(r)" model based on the radius
		// time(r) = c0 + (c1 * r)
		double c0 = 35.9225; //Not very predictable. Simply used average.
		double c1 = tprCoefficient;

		
		// Formula for the area under "time(r)" for range r_start to r_end (trapezoid area)
		// time(r1, r2) = (1/2) * (range) * (r1 + r2)
		double timeStart = c0 + (c1 * r_start);
		double timeEnd = c0 + (c1 * r_end);
		double timeEstimate = (0.5) * (r_end - r_start) * (timeStart + timeEnd);

		// Improve fit for large values of N
		// accuracy within 1% of actual values at N > 300
		double d0 = -0.16741758251155;
		double d1 = 0.209443008293443;
		double correctionCoefficient = d0 + d1 * Math.log(n);
		
		//Apply correction
		double correctedTimeEstimate = correctionCoefficient * timeEstimate;
		
		return correctedTimeEstimate;
	}
	
	/**
	 * Draws circles defined by a list of maxima on the input image
	 * @param draw
	 * @param maxima
	 * @param w
	 * @param h
	 * @param file
	 * @throws IOException
	 */
	private void drawOverlay(BufferedImage draw, List<Integer[]> maxima, int w, int h, File file) throws IOException {
		for(Integer[] m : maxima) {
			for (Integer[] p : circleTemplate(m[0])) {
				int rx = m[1] + p[0];
				int ry = m[2] + p[1];
				if (rx < 0) continue;
				if (rx >= w) continue;
				if (ry < 0) continue;
				if (ry >= h) continue;
				draw.setRGB(rx, ry, Color.red.getRGB());
			}
		}
		for(Integer[] m : maxima) {
			draw.setRGB(m[1], m[2], Color.green.getRGB());
		}
		draw.flush();
		ImageIO.write(draw, "png", file);
	}
	
	/**
	 * Memoized generation of circle templates
	 * @param r
	 * @return
	 */
	public List<Integer[]> circleTemplate (int r) {
		//Precalculate circle kernels
		if (!circleKernel.containsKey(r)) {
			circleKernel.put(r, bresenhamCircle(r));
		}
		return circleKernel.get(r);
	}
	
	/**
	 * Draw a circle by drawing points at finely discretized values of theta
	 * @param r
	 * @return
	 */
	public List<Integer[]> objectivelyWorseCircle(int r) {
		double thetaStep = Math.PI / (Math.pow(r, 2));
		List<Integer[]> kernel = new ArrayList<Integer[]>();
		int lastx = -1;
		int lasty = -1;
		for (double theta = 0.0; theta < 2 * Math.PI; theta += thetaStep) {
			int x = (int) ((int)((Math.cos(theta) * (double)r)));
			int y = (int) ((int)((Math.sin(theta) * (double)r)));
			if (x == lastx && y == lasty) continue;
			lastx = x;
			lasty = y;
			kernel.add(point(x, y));
			
		}
		return kernel;
	}
	
	/**
	 * A fast algorithm for drawing circles, courtesy of Wikipedia
	 * https://en.wikipedia.org/wiki/Midpoint_circle_algorithm
	 * @param r
	 * @return
	 */
	public List<Integer[]> bresenhamCircle(int r) {
		int x = r;
		int y = 0;
		int decisionOver2 = 1 - x;
		
		List<Integer[]> kernel = new ArrayList<Integer[]>();
		while (y <= x) {
			kernel.add(point(x, y)); // Octant 1
			kernel.add(point(y, x)); // Octant 2
			kernel.add(point(-x, y)); // Octant 4
			kernel.add(point(-y, x)); // Octant 3
			kernel.add(point(-x, -y)); // Octant 5
			kernel.add(point(-y, -x)); // Octant 6
			kernel.add(point(x, -y)); // Octant 7
			kernel.add(point(y, -x)); // Octant 8
			y++;
			if (decisionOver2 <= 0) {
				decisionOver2 += 2 * y + 1;
			} else {
				x--;
				decisionOver2 += 2 * (y - x) + 1;
			}
		}
		return kernel;
	}
	
	/**
	 * Convert an image to grayscale
	 * @param img
	 * @return
	 */
	private static int[] grayscale(BufferedImage img) {
		byte[] src = ((DataBufferByte) img.getRaster().getDataBuffer()).getData();
		
		int[] gray = new int[src.length/3];
		for (int i = 0; i < gray.length; i ++) {
			/* 
			 * Java byte primitive is signed, whereas the values
			 * we're reading are unsigned. We will store the image
			 * as an int[] array to preserve positive values above
			 * 127.
			 */
			gray[i] = ((src[i*3] & 0xFF) + (src[i*3 + 1] & 0xFF) + (src[i*3 + 2] & 0xFF)) / 3;
		}
		
		return gray;
	}
	public ArrayList<Integer[]> hough(BufferedImage img, int r_start, int r_end) throws Exception {
		int w = img.getWidth();
		int h = img.getHeight();
		int r_total = r_end - r_start;
		if (r_start < 5) r_start = 5;
		// Create boolean mask to work with
		int[] gray = grayscale(img);
		
		// Basic thresholding
		boolean[] thresh = new boolean[gray.length];
		for (int i = 0; i < gray.length; i++) {
			thresh[i] = gray[i] < 128;
		}
		// Increment accumulator for all circles for all edge points
		int[][][] accumulator = new int[r_total][w][h];
		for (int r = r_start; r < r_end; r++) {
			//long time = System.currentTimeMillis();
			for (int y = 0; y < h; y++) {
				for (int x = 0; x < w; x++) {
					// For all edge pixels
					if (thresh[(y * w) + x]) {
						for (Integer[] p : circleTemplate(r)) {
							int rx = x + p[0];
							int ry = y + p[1];
							if (rx < 0) continue;
							if (rx >= w) continue;
							if (ry < 0) continue;
							if (ry >= h) continue;
							accumulator[r - r_start][rx][ry] ++;
						}
					}
				}
			}
			//System.out.println(w + "\t" + r + "\t" + (System.currentTimeMillis() - time));
		}

		// Find local maxima
		int threshRadius = 2;
		List<Integer[]> maxima = new ArrayList<Integer[]>();
		for (int r = r_start; r < r_end; r++) {
			for (int x = 0; x < w; x++) {
				for (int y = 0; y < h; y++) {
					
					// Define start and end locations for a box
					// Crop at edges
					int i_start = r - threshRadius;
					int i_end = r + threshRadius;
					if (i_start < r_start) i_start = r_start;
					if (i_end > r_end) i_end = r_end;
					
					int j_start = x - threshRadius;
					int j_end = x + threshRadius;
					if (j_start < 0) j_start = 0;
					if (j_end > w) j_end = w;

					int k_start = y - threshRadius;
					int k_end = y + threshRadius;
					if (k_start < 0) k_start = 0;
					if (k_end > h) k_end = h;
					
					double total = (i_end - i_start) * (j_end - j_start) * (k_end - k_start);
					
					// Find the maximum and sum of elements in this box
					int sum = 0;
					int max = 0;
					for (int i = i_start; i < i_end; i++) {
						for (int j = j_start; j < j_end; j++) {
							for (int k = k_start; k < k_end; k++) {
								sum += accumulator[i - r_start][j][k];
								if (accumulator[i - r_start][j][k] > max) max = accumulator[i - r_start][j][k];
							}
						}
					}
					
					// Decide if it's probably a local maximum
					double avg = sum / total;
					int val = accumulator[r - r_start][x][y];
					if ( val == max && val > avg) {
						maxima.add(maximum(r, x, y, val));
					}
				}
			}
		}
		
		// Prune unlikely maxima
		ArrayList<Integer[]> likelyMaxima = new ArrayList<Integer[]>();
		for (Integer[] m : maxima) {
			double fullCircleValue = circleTemplate(m[0]).size();
			// Check that this maximum contains at least half a circle's worth of points...
			if (m[3] / fullCircleValue > 0.9) {
				likelyMaxima.add(m);
			}
		}
		
		return likelyMaxima;
		//drawOverlay(thresh, likelyMaxima, w, h, "overlay.png");
	}
	
	private static Integer[] point(int x, int y) {
		return new Integer[] {x, y};
	}

	private static Integer[] maximum(int r, int x, int y, int v) {
		return new Integer[] {r, x, y, v};
	}
}
