//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Implementation of the HEALPix and rHEALPix projections.
//			For background see <http://code.scenzgrid.org/index.php/p/scenzgrid-py/source/tree/master/docs/rhealpix_dggs.pdf>.
// Authors:	Alex Raichev (raichev@cs.auckland.ac.nz) 
//			Michael Speth (spethm@landcareresearch.co.nz)
// Notes:	Raichev implemented these projections in Python and 
//			Speth translated them into C here.
//
//*****************************************************************************
// Copyright (c) 2001, Thomas Flemming, tf@ttqv.com
// Copyright (c) 2011-2015 by the Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//****************************************************************************

using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_healpix : PJ
	{
		protected double qp;
		protected double[] apa;

		public override string Name { get { return "healpix"; } }
		public override string DescriptionName { get { return "HEALPix"; } }
		public override string DescriptionType { get { return "Sph&Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		static double[][] R1=new double[][] { new double[] { 0, -1 }, new double[] { 1, 0 } }; // Matrix for counterclockwise rotation by pi/2
		static double[][] R2=new double[][] { new double[] { -1, 0 }, new double[] { 0, -1 } }; // Matrix for counterclockwise rotation by pi
		static double[][] R3=new double[][] { new double[] { 0, 1 }, new double[] { -1, 0 } }; // Matrix for counterclockwise rotation by 3*pi/2
		static double[][] IDENT=new double[][] { new double[] { 1, 0 }, new double[] { 0, 1 } }; // Identity matrix

		// IDENT, R1, R2, R3, R1 inverse, R2 inverse, R3 inverse:
		static double[][][] ROT=new double[][][] { IDENT, R1, R2, R3, R3, R2, R1 };

		// Fuzz to handle rounding errors:
		const double EPS=1e-15;

		enum Region { north, south, equatorial }

		struct CapMap
		{
			public int cn; // An integer 0--3 indicating the position of the polar cap.
			public double x, y; // Coordinates of the pole point (point of most extreme latitude on the polar caps).
			public Region region;
		}

		static double[][][] rot=ROT;

		// Returns the sign of the double.
		// @param v the parameter whose sign is returned.
		// @return 1 for positive number, -1 for negative, and 0 for zero.
		static double pj_sign(double v)
		{
			return v>0?1:(v<0?-1:0);
		}

		// Return the index of the matrix in ROT.
		// @param index ranges from -3 to 3.
		static int get_rotate_index(int index)
		{
			switch(index)
			{
				case 0: return 0;
				case 1: return 1;
				case 2: return 2;
				case 3: return 3;
				case -1: return 4;
				case -2: return 5;
				case -3: return 6;
			}
			return 0;
		}

		// Return 1 if point (testx, testy) lies in the interior of the polygon
		// determined by the vertices in vert, and return 0 otherwise.
		// See http://paulbourke.net/geometry/polygonmesh/ for more details.
		// @param nvert the number of vertices in the polygon.
		// @param vert the (x, y)-coordinates of the polygon's vertices
		static bool pnpoly(int nvert, double[][] vert, double testx, double testy)
		{
			// Check for boundrary cases
			for(int i=0; i<nvert; i++)
			{
				if(testx==vert[i][0]&&testy==vert[i][1]) return true;
			}

			double x1=vert[0][0];
			double y1=vert[0][1];
			double x2, y2;

			int counter=0;
			for(int i=1; i<nvert; i++, x1=x2, y1=x2)
			{
				x2=vert[i%nvert][0];
				y2=vert[i%nvert][1];

				if(testy<=Math.Min(y1, y2)||testy>Math.Max(y1, y2)||testx>Math.Max(x1, x2)) continue;

				if(y1==y2) continue;

				double xinters=(testy-y1)*(x2-x1)/(y2-y1)+x1;
				if(x1==x2||testx<=xinters) counter++;
			}

			return counter%2!=0;
		}

		// Return 1 if (x, y) lies in (the interior or boundary of) the image of the
		// HEALPix projection (in case proj=0) or in the image the rHEALPix projection
		// (in case proj=1), and return 0 otherwise.
		// @param north_square the position of the north polar square (rHEALPix only)
		// @param south_square the position of the south polar square (rHEALPix only)
		protected static bool in_image(double x, double y, int proj, int north_square, int south_square)
		{
			if(proj==0)
			{
				double[][] healpixVertsJit=new double[][]
				{
					new double[] { -1.0*Math.PI-EPS, Math.PI/4.0 },
					new double[] { -3.0*Math.PI/4.0, Math.PI/2.0+EPS },
					new double[] { -1.0*Math.PI/2.0, Math.PI/4.0+EPS },
					new double[] { -1.0*Math.PI/4.0, Math.PI/2.0+EPS },
					new double[] { 0.0, Math.PI/4.0+EPS },
					new double[] { Math.PI/4.0, Math.PI/2.0+EPS },
					new double[] { Math.PI/2.0, Math.PI/4.0+EPS },
					new double[] { 3.0*Math.PI/4.0, Math.PI/2.0+EPS },
					new double[] { Math.PI+EPS, Math.PI/4.0 },
					new double[] { Math.PI+EPS, -1.0*Math.PI/4.0 },
					new double[] { 3.0*Math.PI/4.0, -1.0*Math.PI/2.0-EPS },
					new double[] { Math.PI/2.0, -1.0*Math.PI/4.0-EPS },
					new double[] { Math.PI/4.0, -1.0*Math.PI/2.0-EPS },
					new double[] { 0.0, -1.0*Math.PI/4.0-EPS },
					new double[] { -1.0*Math.PI/4.0, -1.0*Math.PI/2.0-EPS },
					new double[] { -1.0*Math.PI/2.0, -1.0*Math.PI/4.0-EPS },
					new double[] { -3.0*Math.PI/4.0, -1.0*Math.PI/2.0-EPS },
					new double[] { -1.0*Math.PI-EPS, -1.0*Math.PI/4.0 }
				};
				return pnpoly(healpixVertsJit.Length, healpixVertsJit, x, y);
			}
			else
			{
				double[][] rhealpixVertsJit=new double[][]
				{
					new double[] { -1.0*Math.PI-EPS, Math.PI/4.0+EPS },
					new double[] { -1.0*Math.PI+north_square*Math.PI/2.0-EPS, Math.PI/4.0+EPS },
					new double[] { -1.0*Math.PI+north_square*Math.PI/2.0-EPS, 3*Math.PI/4.0+EPS },
					new double[] { -1.0*Math.PI+(north_square+1.0)*Math.PI/2.0+EPS, 3*Math.PI/4.0+EPS },
					new double[] { -1.0*Math.PI+(north_square+1.0)*Math.PI/2.0+EPS, Math.PI/4.0+EPS },
					new double[] { Math.PI+EPS, Math.PI/4.0+EPS },
					new double[] { Math.PI+EPS, -1.0*Math.PI/4.0-EPS },
					new double[] { -1.0*Math.PI+(south_square+1.0)*Math.PI/2.0+EPS, -1.0*Math.PI/4.0-EPS },
					new double[] { -1.0*Math.PI+(south_square+1.0)*Math.PI/2.0+EPS, -3.0*Math.PI/4.0-EPS },
					new double[] { -1.0*Math.PI+south_square*Math.PI/2.0-EPS, -3.0*Math.PI/4.0-EPS },
					new double[] { -1.0*Math.PI+south_square*Math.PI/2.0-EPS, -1.0*Math.PI/4.0-EPS },
					new double[] { -1.0*Math.PI-EPS, -1.0*Math.PI/4.0-EPS }
				};
				return pnpoly(rhealpixVertsJit.Length, rhealpixVertsJit, x, y);
			}
		}

		// Return the authalic latitude of latitude alpha (if inverse=0) or
		// return the approximate latitude of authalic latitude alpha (if inverse=1).
		// P contains the relavent ellipsoid parameters.
		protected static double auth_lat(PJ_healpix P, double alpha, bool inverse)
		{
			if(!inverse)
			{
				// Authalic latitude.
				double q=Proj.pj_qsfn(Math.Sin(alpha), P.e, 1.0-P.es);
				double qp=P.qp;
				double ratio=q/qp;
				if(Math.Abs(ratio)>1) ratio=pj_sign(ratio); // Rounding error.
				return Math.Asin(ratio);
			}

			// Approximation to inverse authalic latitude.
			return Proj.pj_authlat(alpha, P.apa);
		}

		// Return the HEALPix projection of the longitude-latitude point lp on
		// the unit sphere.
		protected XY healpix_sphere(LP lp)
		{
			double lam=lp.lam;
			double phi=lp.phi;
			double phi0=Math.Asin(2.0/3.0);

			XY xy;

			// equatorial region
			if(Math.Abs(phi)<=phi0)
			{
				xy.x=lam;
				xy.y=3.0*Math.PI/8.0*Math.Sin(phi);
			}
			else
			{
				double sigma=Math.Sqrt(3.0*(1-Math.Abs(Math.Sin(phi))));
				double cn=Math.Floor(2*lam/Math.PI+2);
				if(cn>=4) cn=3;

				double lamc=-3*Math.PI/4+(Math.PI/2)*cn;
				xy.x=lamc+(lam-lamc)*sigma;
				xy.y=pj_sign(phi)*Math.PI/4*(2-sigma);
			}

			return xy;
		}

		// Return the inverse of healpix_sphere().
		protected LP healpix_sphere_inverse(XY xy)
		{
			double x=xy.x;
			double y=xy.y;
			double y0=Math.PI/4.0;

			LP lp;

			// Equatorial region.
			if(Math.Abs(y)<=y0)
			{
				lp.lam=x;
				lp.phi=Math.Asin(8.0*y/(3.0*Math.PI));
			}
			else if(Math.Abs(y)<Math.PI/2.0)
			{
				double cn=Math.Floor(2.0*x/Math.PI+2.0);
				if(cn>=4) cn=3;

				double xc=-3.0*Math.PI/4.0+(Math.PI/2.0)*cn;
				double tau=2.0-4.0*Math.Abs(y)/Math.PI;
				lp.lam=xc+(x-xc)/tau;
				lp.phi=pj_sign(y)*Math.Asin(1.0-Math.Pow(tau, 2.0)/3.0);
			}
			else
			{
				lp.lam=-1.0*Math.PI;
				lp.phi=pj_sign(y)*Math.PI/2.0;
			}

			return lp;
		}

		// Return the vector sum a + b, where a and b are 2-dimensional vectors.
		// @param ret holds a + b.
		static void vector_add(double[] a, double[] b, double[] ret)
		{
			for(int i=0; i<2; i++) ret[i]=a[i]+b[i];
		}

		// Return the vector difference a - b, where a and b are 2-dimensional vectors.
		// @param ret holds a - b.
		static void vector_sub(double[] a, double[] b, double[] ret)
		{
			for(int i=0; i<2; i++) ret[i]=a[i]-b[i];
		}

		// Return the 2 x 1 matrix product a*b, where a is a 2 x 2 matrix and 
		// b is a 2 x 1 matrix.
		// @param ret holds a*b.
		static void dot_product(double[][] a, double[] b, double[] ret)
		{
			int length=2;
			for(int i=0; i<length; i++)
			{
				ret[i]=0;
				for(int j=0; j<length; j++) ret[i]+=a[i][j]*b[j];
			}
		}

		// Return the number of the polar cap, the pole point coordinates, and the region that (x, y) lies in.
		// If inverse=0, then assume (x,y) lies in the image of the HEALPix projection of the unit sphere.
		// If inverse=1, then assume (x,y) lies in the image of the (north_square, south_square)-rHEALPix projection of the unit sphere.
		static CapMap get_cap(double x, double y, int north_square, int south_square, bool inverse)
		{
			CapMap capmap;
			capmap.x=x;
			capmap.y=y;

			if(!inverse)
			{
				double c;
				if(y>Math.PI/4.0)
				{
					capmap.region=Region.north;
					c=Math.PI/2.0;
				}
				else if(y<-1*Math.PI/4.0)
				{
					capmap.region=Region.south;
					c=-1*Math.PI/2.0;
				}
				else
				{
					capmap.region=Region.equatorial;
					capmap.cn=0;
					return capmap;
				}

				capmap.y=c;

				// polar region
				if(x<-1*Math.PI/2.0)
				{
					capmap.cn=0;
					capmap.x=(-1*3.0*Math.PI/4.0);
				}
				else if(x>=-1*Math.PI/2.0&&x<0)
				{
					capmap.cn=1;
					capmap.x=-1*Math.PI/4.0;
				}
				else if(x>=0&&x<Math.PI/2.0)
				{
					capmap.cn=2;
					capmap.x=Math.PI/4.0;
				}
				else
				{
					capmap.cn=3;
					capmap.x=3.0*Math.PI/4.0;
				}

				return capmap;
			}
			else
			{
				// Polar Region, find the HEALPix polar cap number that
				// x, y moves to when rHEALPix polar square is disassembled.
				double eps=1e-15; // Kludge. Fuzz to avoid some rounding errors.

				if(y>Math.PI/4.0)
				{
					capmap.region=Region.north;
					capmap.x=-3.0*Math.PI/4.0+north_square*Math.PI/2.0;
					capmap.y=Math.PI/2.0;
					x=x-north_square*Math.PI/2.0;

					if(y>=-1*x-Math.PI/4.0-eps&&y<x+5.0*Math.PI/4.0-eps) capmap.cn=(north_square+1)%4;
					else if(y>-1*x-1*Math.PI/4.0+eps&&y>=x+5.0*Math.PI/4.0-eps) capmap.cn=(north_square+2)%4;
					else if(y<=-1*x-1*Math.PI/4.0+eps&&y>x+5.0*Math.PI/4.0+eps) capmap.cn=(north_square+3)%4;
					else capmap.cn=north_square;
				}
				else if(y<-1*Math.PI/4.0)
				{
					capmap.region=Region.south;
					capmap.x=-3.0*Math.PI/4.0+south_square*Math.PI/2;
					capmap.y=-1*Math.PI/2.0;
					x=x-south_square*Math.PI/2.0;

					if(y<=x+Math.PI/4.0+eps&&y>-1*x-5.0*Math.PI/4.0+eps) capmap.cn=(south_square+1)%4;
					else if(y<x+Math.PI/4.0-eps&&y<=-1*x-5.0*Math.PI/4.0+eps) capmap.cn=(south_square+2)%4;
					else if(y>=x+Math.PI/4.0-eps&&y<-1*x-5.0*Math.PI/4.0-eps) capmap.cn=(south_square+3)%4;
					else capmap.cn=south_square;
				}
				else
				{
					capmap.region=Region.equatorial;
					capmap.cn=0;
				}

				return capmap;
			}
		}

		// Rearrange point (x, y) in the HEALPix projection by combining the polar caps into two polar squares.
		// Put the north polar square in position north_square and the south polar square in position south_square.
		// If inverse=1, then uncombine the polar caps.
		// @param north_square integer between 0 and 3.
		// @param south_square integer between 0 and 3.
		protected static XY combine_caps(double x, double y, int north_square, int south_square, bool inverse)
		{
			CapMap capmap=get_cap(x, y, north_square, south_square, inverse);

			XY xy;
			if(capmap.region==Region.equatorial)
			{
				xy.x=capmap.x;
				xy.y=capmap.y;
				return xy;
			}

			double[][] tmpRot;
			double[] a=new double[2];

			if(!inverse)
			{
				// Rotate (x, y) about its polar cap tip and then translate it to
				// north_square or south_square.

				if(capmap.region==Region.north)
				{
					int pole=north_square;
					a[0]=-3.0*Math.PI/4.0+pole*Math.PI/2;
					a[1]=Math.PI/2.0+pole*0;
					tmpRot=rot[get_rotate_index(capmap.cn-pole)];
				}
				else
				{
					int pole=south_square;
					a[0]=-3.0*Math.PI/4.0+pole*Math.PI/2;
					a[1]=Math.PI/-2.0+pole*0;
					tmpRot=rot[get_rotate_index(-1*(capmap.cn-pole))];
				}
			}
			else
			{
				// Inverse function.
				// Unrotate (x, y) and then translate it back.

				// disassemble
				if(capmap.region==Region.north)
				{
					int pole=north_square;
					a[0]=-3.0*Math.PI/4.0+capmap.cn*Math.PI/2;
					a[1]=Math.PI/2.0+capmap.cn*0;
					tmpRot=rot[get_rotate_index(-1*(capmap.cn-pole))];
				}
				else
				{
					int pole=south_square;
					a[0]=-3.0*Math.PI/4.0+capmap.cn*Math.PI/2;
					a[1]=Math.PI/-2.0+capmap.cn*0;
					tmpRot=rot[get_rotate_index(capmap.cn-pole)];
				}
			}

			// translate, rotate, then shift
			double[] v=new double[2] { x, y };
			double[] c=new double[2] { capmap.x, capmap.y };

			double[] v_min_c=new double[2];
			vector_sub(v, c, v_min_c);

			double[] ret_dot=new double[2];
			dot_product(tmpRot, v_min_c, ret_dot);

			double[] vector=new double[2];
			vector_add(ret_dot, a, vector);

			xy.x=vector[0];
			xy.y=vector[1];

			return xy;
		}

		// spheroid
		XY s_healpix_forward(LP lp)
		{
			return healpix_sphere(lp);
		}

		// ellipsoid
		XY e_healpix_forward(LP lp)
		{
			lp.phi=auth_lat(this, lp.phi, false);
			return healpix_sphere(lp);
		}

		// spheroid
		LP s_healpix_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			// Check whether (x, y) lies in the HEALPix image
			if(!in_image(xy.x, xy.y, 0, 0, 0))
			{
				lp.lam=Libc.HUGE_VAL;
				lp.phi=Libc.HUGE_VAL;
				Proj.pj_ctx_set_errno(ctx, -15);

				return lp;
			}

			return healpix_sphere_inverse(xy);
		}

		// ellipsoid
		LP e_healpix_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			// Check whether (x, y) lies in the HEALPix image.
			if(!in_image(xy.x, xy.y, 0, 0, 0))
			{
				lp.lam=Libc.HUGE_VAL;
				lp.phi=Libc.HUGE_VAL;
				Proj.pj_ctx_set_errno(ctx, -15);

				return lp;
			}

			lp=healpix_sphere_inverse(xy);

			lp.phi=auth_lat(this, lp.phi, true);

			return lp;
		}

		public override PJ Init()
		{
			if(es!=0)
			{
				apa=Proj.pj_authset(es); // For auth_lat().
				qp=Proj.pj_qsfn(1.0, e, one_es); // For auth_lat().
				a=a*Math.Sqrt(0.5*qp); // Set a to authalic radius.
				ra=1.0/a;
				fwd=e_healpix_forward;
				inv=e_healpix_inverse;
			}
			else
			{
				fwd=s_healpix_forward;
				inv=s_healpix_inverse;
			}

			return this;
		}
	}

	class PJ_rhealpix : PJ_healpix
	{
		protected int north_square;
		protected int south_square;

		public override string Name { get { return "rhealpix"; } }
		public override string DescriptionName { get { return "rHEALPix"; } }
		public override string DescriptionType { get { return "Sph&Ell"; } }
		public override string DescriptionParameters { get { return "north_square= south_square="; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(" +north_square={0}", north_square);
				ret.AppendFormat(" +south_square={0}", south_square);
				return ret.ToString();
			}
		}

		// spheroid
		XY s_rhealpix_forward(LP lp)
		{
			XY xy=healpix_sphere(lp);
			return combine_caps(xy.x, xy.y, north_square, south_square, false);
		}

		// ellipsoid
		XY e_rhealpix_forward(LP lp)
		{
			lp.phi=auth_lat(this, lp.phi, false);
			XY xy=healpix_sphere(lp);
			return combine_caps(xy.x, xy.y, north_square, south_square, false);
		}

		// spheroid
		LP s_rhealpix_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			// Check whether (x, y) lies in the rHEALPix image.
			if(!in_image(xy.x, xy.y, 1, north_square, south_square))
			{
				lp.lam=Libc.HUGE_VAL;
				lp.phi=Libc.HUGE_VAL;
				Proj.pj_ctx_set_errno(ctx, -15);
				return lp;
			}

			xy=combine_caps(xy.x, xy.y, north_square, south_square, true);

			return healpix_sphere_inverse(xy);
		}

		// ellipsoid
		LP e_rhealpix_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			// Check whether (x, y) lies in the rHEALPix image.
			if(!in_image(xy.x, xy.y, 1, north_square, south_square))
			{
				lp.lam=Libc.HUGE_VAL;
				lp.phi=Libc.HUGE_VAL;
				Proj.pj_ctx_set_errno(ctx, -15);
				return lp;
			}

			xy=combine_caps(xy.x, xy.y, north_square, south_square, true);
			lp=healpix_sphere_inverse(xy);
			lp.phi=auth_lat(this, lp.phi, true);

			return lp;
		}

		public override PJ Init()
		{
			north_square=Proj.pj_param_i(ctx, parameters, "north_square");
			south_square=Proj.pj_param_i(ctx, parameters, "south_square");

			// Check for valid north_square and south_square inputs.
			if(north_square<0||north_square>3)
			{
				Proj.pj_ctx_set_errno(ctx, -47);
				return null;
			}

			if(south_square<0||south_square>3)
			{
				Proj.pj_ctx_set_errno(ctx, -47);
				return null;
			}

			if(es!=0)
			{
				apa=Proj.pj_authset(es); // For auth_lat().
				qp=Proj.pj_qsfn(1.0, e, one_es); // For auth_lat().
				a=a*Math.Sqrt(0.5*qp); // Set a to authalic radius.
				ra=1.0/a;
				fwd=e_rhealpix_forward;
				inv=e_rhealpix_inverse;
			}
			else
			{
				fwd=s_rhealpix_forward;
				inv=s_rhealpix_inverse;
			}

			return this;
		}
	}
}
