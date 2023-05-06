package triangulation;

import static java.lang.Math.cos;
import static java.lang.Math.sin;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.StringTokenizer;

import fem.Element;
import fem.Model;
import fem.Node;
import math.Vect;
import math.util;

class Point implements Comparable<Point>{
	double x, y;

	public int compareTo(Point p) {
	double r=0;
		if (this.x == p.x) {
			r=this.y- p.y;

		} else {
			r=this.x- p.x;
		}
		if(r>0) return 1;
		else if (r==0) return 0;
		else return -1;
	}

	public String toString() {
		return "(" + x + "," + y + ")\n";
	}

}

public class ConvexHull {

	public static double cross(Point O, Point A, Point B) {
		return (A.x - O.x) * (double) (B.y - O.y) - (A.y - O.y) * (double) (B.x - O.x);
	}

	public static Point[] convex_hull(Point[] P) {

		if (P.length > 1) {
			int n = P.length, k = 0;
			Point[] H = new Point[2 * n];

			Arrays.sort(P);

			// Build lower hull
			for (int i = 0; i < n; ++i) {
				while (k >= 2 && cross(H[k - 2], H[k - 1], P[i]) <= 0)
					k--;
				H[k++] = P[i];
			}

			// Build upper hull
			for (int i = n - 2, t = k + 1; i >= 0; i--) {
				while (k >= t && cross(H[k - 2], H[k - 1], P[i]) <= 0)
					k--;
				H[k++] = P[i];
			}
			if (k > 1) {
				H = Arrays.copyOfRange(H, 0, k - 1); // remove non-hull vertices after k; remove k - 1 which is a duplicate
			}
			return H;
		} else if (P.length <= 1) {
			return P;
		} else {
			return null;
		}
	}
	
	public static Vect[] getConvexHull(Vect[] p){
	
		Point[] set = new Point[p.length];
		
		for(int i=0;i<p.length;i++){
			set[i]=new Point();
			set[i].x=p[i].el[0];
			set[i].y=p[i].el[1];
		}
	
		Point[] hull0= convex_hull(set).clone();
		Vect[] hull=new Vect[hull0.length];
		for (int i = 0; i < hull.length; i++) {
			hull[i]=new Vect(hull0[i].x,hull0[i].y);
			
		}
		return hull;
	}

	public static void main(String[] args) throws IOException {

		String bun=util.getFile(0);

		Model model=new Model(bun);
		
		String folder=new File(model.meshFilePath).getParentFile().getPath();

		int nNodes=model.numberOfNodes;
		Vect[] p = new Vect[nNodes];
		
		for(int i=1;i<=nNodes;i++){
			p[i-1]=model.node[i].getCoord();
			//p[i-1].x=Math.random();
			//p[i-1].y=Math.random();
			//p[i-1]=model.node[i].getCoord();
			//p[i-1].x=model.node[i].getCoord(0);//(int)(1000*i*.001*cos(.1*i));
			//p[i-1].y=model.node[i].getCoord(1);
			//P[i-1].el[1]=Math.random();
			//P[i-1].hshow();
			//util.pr(p[i-1].x+"  "+p[i-1].y);
		}
	
		Vect[] hull =getConvexHull(p);

		for (int i = 0; i < hull.length; i++) {
		//	if (hull[i] != null)
			//	System.out.print(hull[i]);
			
		}
		
		Vect center=new Vect(2);
		for (int i = 0; i < hull.length; i++) {
			center=center.add(hull[i]);
		}
		
		center.timesVoid(1./hull.length);

		Vect[] hull2 = new Vect[hull.length];
		for (int i = 0; i < hull.length; i++) {
			hull2[i]=hull[i].sub(center).times(1.01).add(center);

			//System.out.print(hull2[i]);
		}

		int nRegions=1;
		Model hullMesh=new Model(nRegions,hull.length-1,2*hull.length,"quadrangle");

		for(int i=1;i<=nRegions;i++){


			hullMesh.region[i].setMaterial("mat"+i);
			hullMesh.region[i].setName("reg"+i);

		}
	//	for(int i=1;i<=hullMesh.numberOfElements;i++)
		for(int i=0;i<hull.length;i++)
			hullMesh.node[i+1].setCoord(hull[i]);
		for(int i=hull.length;i<2*hull.length;i++)
			hullMesh.node[i+1].setCoord(hull2[i-hull.length]);

		for(int i=1;i<=hullMesh.numberOfElements;i++){
			int[] vn=new int[4];
			if(i<hullMesh.numberOfElements){
			vn[0]=i;
			vn[1]=i+1;
			vn[2]=hull.length+i+1;
			vn[3]=hull.length+i;
			}else{
				vn[0]=i;
				vn[1]=1;
				vn[2]=hull.length+1;
				vn[3]=hull.length+i;	
			}
			hullMesh.element[i].setVertNumb(vn);
		}
		
		hullMesh.region[1].setFirstEl(1);
		hullMesh.region[1].setLastEl(hullMesh.numberOfElements);

		
		hullMesh.writeMesh(folder+"\\hull.txt");
	}

}