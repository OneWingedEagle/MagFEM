package triangulation;
import fem.Edge;
import fem.Node;

class Triangle
{

	Node p1, p2, p3;
	Edge e1,e2,e3;
	boolean isBad;

	public static double epsilon=1e-12;

	Triangle(Node p1, Node p2, Node p3)
		{
			this.p1=p1;
			this.p2=p2;
			this.p3=p3;
			this.e1=new Edge(p1,p2);
			this.e2=new Edge(p2,p3);
			this.e3=new Edge(p3,p1);
			isBad=false;

		}

		boolean containsVertex(final Node v) 
		{
			return almost_equal(p1, v) || almost_equal(p2, v) || almost_equal(p3, v);
		}

		public boolean circumCircleContains(final Node v) 
		{
			final double ab = p1.getCoord().norm2();
			final double cd = p2.getCoord().norm2();
			final double ef = p3.getCoord().norm2();

			final double circum_x = (ab * (p3.getCoord(1) - p2.getCoord(1)) + cd * (p1.getCoord(1) - p3.getCoord(1)) + ef * (p2.getCoord(1) - p1.getCoord(1))) / (p1.getCoord(0) * (p3.getCoord(1) - p2.getCoord(1)) + p2.getCoord(0) * (p1.getCoord(1) - p3.getCoord(1)) + p3.getCoord(0) * (p2.getCoord(1) - p1.getCoord(1)));
			final double circum_y = (ab * (p3.getCoord(0) - p2.getCoord(0)) + cd * (p1.getCoord(0) - p3.getCoord(0)) + ef * (p2.getCoord(0) - p1.getCoord(0))) / (p1.getCoord(1) * (p3.getCoord(0) - p2.getCoord(0)) + p2.getCoord(1) * (p1.getCoord(0) - p3.getCoord(0)) + p3.getCoord(1) * (p2.getCoord(0) - p1.getCoord(0)));

			Node circum=new Node(-1,2);
			circum.setCoord(0,half(circum_x));
			circum.setCoord(1,half(circum_y));
		
			final double circum_radius = p1.dist2(circum);
			final double dist2 = v.dist2(circum);
			return dist2 <= circum_radius;
		}




public static boolean almost_equal(final Triangle t1, final Triangle t2)
{
	return	(almost_equal(t1.p1 , t2.p1) || almost_equal(t1.p1 , t2.p2) || almost_equal(t1.p1 , t2.p3)) &&
			(almost_equal(t1.p2 , t2.p1) || almost_equal(t1.p2 , t2.p2) || almost_equal(t1.p2 , t2.p3)) &&
			(almost_equal(t1.p3 , t2.p1) || almost_equal(t1.p3 , t2.p2) || almost_equal(t1.p3 , t2.p3));
}

public static boolean almost_equal (final Edge e1, final Edge e2)
{
	return	(almost_equal(e1.node[0], e2.node[0]) && almost_equal(e1.node[1], e2.node[1])) ||
			(almost_equal(e1.node[0], e2.node[1]) && almost_equal(e1.node[1], e2.node[0]));
}

public static boolean almost_equal(final Node n1, final Node n2)
{
	return	(n1.getCoord().sub(n2.getCoord()).norm2()<epsilon);
}


private double half(double x){return 0.5 * x;}

}
