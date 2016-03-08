package org.kylemoy.Hough;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.DataBufferByte;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import javax.imageio.ImageIO;

import psimj.Communicator;
import psimj.Task;

public class HoughExceedinglyBroken implements Task {
	
	static String path = "D:/data/hough/circles1024.jpg";
	
	public static void main(String[] args) throws Exception {
		//Pool pool = new Pool("127.0.0.1", 8191);
		//pool.useAuthentication(new PoolKey("./pseudo.key"));
		HoughExceedinglyBroken h = new HoughExceedinglyBroken();
		System.out.println(h.circleTemplate(1024).size());
		/*
		BufferedImage img = ImageIO.read(new File(path));
		Hough h = new Hough();
		h.hough(img, 0, 5);
		*/
		/*
		Communicator comm = new LocalCommunicator(12, new Topology.Switch());
		//Communicator comm = pool.requestCommunicator(-1);
		if (comm != null) {
			comm.runTask(Hough.class);
		}
		*/
		/*
				
		//PrintWriter out = new PrintWriter(new File("circlegrowth.txt"));
		
		
		for (int i = 16; i < 1024; i *= 2) {
			System.out.println("For n = " + i);
			int r_max = (int)(Math.ceil(Math.sqrt(Math.pow(i, 2) + Math.pow(i, 2)))) / 2;
			int r_start = 5;
			int r_end = r_max;
			int n = i;
			int nprocs = 4;
			
			//
			long time = System.currentTimeMillis();
			h.hough(new BufferedImage(i,i,BufferedImage.TYPE_INT_RGB),r_start,r_end);
			int estimate = (int)h.estimateCost(i, r_start, r_end);
			long total = (System.currentTimeMillis() - time);
			System.out.println(i + "\t" + estimate + "\t" + total + "\t" + (((total - estimate) + estimate) / (double)estimate));
			//
			double estimatedTotalCost = h.estimateCost(n, r_start, r_end);
			int costPerNode = (int)Math.ceil(estimatedTotalCost / (float)nprocs);
			int start = r_start;
			System.out.println("Total: " + (r_end - r_start));
			int total = 0;
			for (int j = 0; j < nprocs - 1; j++) {
				int end = (int) h.solveRange(n, costPerNode, start);
				System.out.println("Node " + j + " gets " + start + " to " + end + "(" + (end - start) + ")");
				total += (end - start);
				start = end;
			}
			System.out.println("Node " + (nprocs - 1) + " gets " + start + " to " + r_end + "(" + (r_end - start) + ")");
			total += (r_end - start);
			System.out.println("Cumulative: " + total + "\tComplete?:" + (total == (r_end - r_start)));
		}
		
		//out.close();
		//h.hough(img, 0, 128);
		 */
	}

	HashMap<Integer,List<Point>> circleKernel = new HashMap<Integer, List<Point>>();
	int r_start;
	int r_end;
	
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
			int r_max = (int)(Math.ceil(Math.sqrt(Math.pow(w, 2) + Math.pow(h, 2)))) / 2;
			//r_max *= 0.8;
			
			int rValuesPerNode = (int)Math.ceil((r_max / (double)comm.nprocs()));
			
			int start = comm.rank() * rValuesPerNode;
			int end = (comm.rank() + 1) * rValuesPerNode;
			if (comm.rank() == 0) {
				for (int i = 0; i < comm.nprocs(); i++) {
					System.out.println("Rank " + i + " gets rValues " + (i * rValuesPerNode) + " to " + ((i + 1) * rValuesPerNode));
				}
			}
			ArrayList<Maximum> maxima = hough(img, start, end);
			System.out.println("Node " + comm.rank() + " complete.");
			List<ArrayList<Maximum>> maximaList = (List<ArrayList<Maximum>>) comm.all2one_collect(0, maxima, maxima.getClass());
			
			if (comm.rank() == 0) {
				ArrayList<Maximum> allMaxima = new ArrayList<Maximum>();
				for (ArrayList<Maximum> list : maximaList) {
					for (Maximum m : list) {
						allMaxima.add(m);
					}
				}
				System.out.println("Maxima found: " + allMaxima.size());
				File f = new File("overlay.png");
				drawOverlay(img, allMaxima, w, h, f);
				System.out.println("Overlay saved as " + f.getAbsolutePath());
				System.out.println("Total elapsed time: " + (System.currentTimeMillis() - time) + "ms");
			}
			
			comm.finish();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public List<Point> circleTemplate (int r) {
		//Precalculate circle kernels
		if (!circleKernel.containsKey(r)) {
			double thetaStep = Math.PI / (Math.pow(r, 2));
			List<Point> kernel = new ArrayList<Point>();
			int lastx = -1;
			int lasty = -1;
			for (double theta = 0.0; theta < 2 * Math.PI; theta += thetaStep) {
				int x = (int) ((int)((Math.cos(theta) * (double)r)));
				int y = (int) ((int)((Math.sin(theta) * (double)r)));
				if (x == lastx && y == lasty) continue;
				lastx = x;
				lasty = y;
				kernel.add(new Point(x, y));
				
			}
			circleKernel.put((int)r, kernel);
		}
		return circleKernel.get(r);
	}
	public ArrayList<Maximum> hough(BufferedImage img, int start, int end) throws Exception {
		int w = img.getWidth();
		int h = img.getHeight();
		int r_max = (int)(Math.ceil(Math.sqrt(Math.pow(w, 2) + Math.pow(h, 2)))) / 2;
		r_start = start;
		r_end = end;
		int r_total = r_end - r_start;
		if (r_start < 5) r_start = 5;
		if (r_end > r_max) r_end = r_max;
		
		// Create boolean mask to work with
		int[] gray = grayscale(img);
		
		// Basic thresholding
		boolean[] thresh = new boolean[gray.length];
		for (int i = 0; i < gray.length; i++) {
			thresh[i] = gray[i] < 128;
		}
		
		// Increment accumulator for all circles for all edge points
		int[][][] accumulator = new int[r_max][w][h];
		for (int y = 0; y < h; y++) {
			for (int x = 0; x < w; x++) {
				// For all edge pixels
				if (thresh[(y * w) + x]) {
					for (int r = r_start; r < r_end; r++) {
						for (Point p : circleTemplate(r)) {
							int rx = x + p.x;
							int ry = y + p.y;
							if (rx < 0) continue;
							if (rx >= w) continue;
							if (ry < 0) continue;
							if (ry >= h) continue;
							accumulator[r][rx][ry] ++;
						}
					}
				}
			}
		}

		// Find local maxima
		int threshRadius = 2;
		List<Maximum> maxima = new ArrayList<Maximum>();
		for (int r = r_start; r < r_end; r++) {
			for (int x = 0; x < w; x++) {
				for (int y = 0; y < h; y++) {
					
					// Define start and end locations for a box
					// Crop at edges
					int i_start = r - threshRadius;
					int i_end = r + threshRadius;
					if (i_start < 0) i_start = 0;
					if (i_end > r_max) i_end = r_max;
					
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
								sum += accumulator[i][j][k];
								if (accumulator[i][j][k] > max) max = accumulator[i][j][k];
							}
						}
					}
					
					// Decide if it's probably a local maximum
					double avg = sum / total;
					int val = accumulator[r][x][y];
					if ( val == max && val > avg) {
						maxima.add(new Maximum(r, x, y, val));
					}
				}
			}
		}
		
		// Prune unlikely maxima
		ArrayList<Maximum> likelyMaxima = new ArrayList<Maximum>();
		for (Maximum m : maxima) {
			double fullCircleValue = circleTemplate(m.r).size();
			// Check that this maximum contains at least half a circle's worth of points...
			if (m.v / fullCircleValue > 0.6) {
				likelyMaxima.add(m);
			}
		}
		
		return likelyMaxima;
		//drawOverlay(thresh, likelyMaxima, w, h, "overlay.png");
	}

	private void drawOverlay(BufferedImage draw, List<Maximum> maxima, int w, int h, File file) throws IOException {
		for(Maximum m : maxima) {
			for (Point p : circleTemplate(m.r)) {
				int rx = m.x + p.x;
				int ry = m.y + p.y;
				if (rx < 0) continue;
				if (rx >= w) continue;
				if (ry < 0) continue;
				if (ry >= h) continue;
				draw.setRGB(rx, ry, Color.red.getRGB());
			}
		}
		for(Maximum m : maxima) {
			draw.setRGB(m.x, m.y, Color.green.getRGB());
		}
		draw.flush();
		ImageIO.write(draw, "png", file);
	}

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
	private static class Point {
		public final int x, y;
		public Point( int _x , int _y ) {
			x = _x;
			y = _y;
		}
	}

	private static class Maximum implements Serializable{
		private static final long serialVersionUID = 5512457969003041131L;
		public final int x, y, r, v;
		public Maximum( int _r, int _x , int _y, int _v) {
			r = _r;
			x = _x;
			y = _y;
			v = _v;
		}
	}
}
