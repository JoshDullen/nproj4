// This implements the Quadrilateralized Spherical Cube (QSC) projection.
//
// Copyright (c) 2011, 2012  Martin Lambers <marlam@marlam.de>
// Copyright (c) 2014-2015 by the Authors
//
// The QSC projection was introduced in:
// [OL76]
// E.M. O'Neill and R.E. Laubscher, "Extended Studies of a Quadrilateralized
// Spherical Cube Earth Data Base", Naval Environmental Prediction Research
// Facility Tech. Report NEPRF 3-76 (CSC), May 1976.
//
// The preceding shift from an ellipsoid to a sphere, which allows to apply
// this projection to ellipsoids as used in the Ellipsoidal Cube Map model,
// is described in
// [LK12]
// M. Lambers and A. Kolb, "Ellipsoidal Cube Maps for Accurate Rendering of
// Planetary-Scale Terrain Data", Proc. Pacfic Graphics (Short Papers), Sep.
// 2012
//
// You have to choose one of the following projection centers,
// corresponding to the centers of the six cube faces:
// phi0 = 0.0, lam0 = 0.0       ("front" face)
// phi0 = 0.0, lam0 = 90.0      ("right" face)
// phi0 = 0.0, lam0 = 180.0     ("back" face)
// phi0 = 0.0, lam0 = -90.0     ("left" face)
// phi0 = 90.0                  ("top" face)
// phi0 = -90.0                 ("bottom" face)
// Other projection centers will not work!
//
// In the projection code below, each cube face is handled differently.
// See the computation of the face parameter in the Init function
// and the handling of different face values (FACE_*) in the forward and
// inverse projections.
//
// Furthermore, the projection is originally only defined for theta angles
// between (-1/4 * PI) and (+1/4 * PI) on the current cube face. This area
// of definition is named AREA_0 in the projection code below. The other
// three areas of a cube face are handled by rotation of AREA_0.

using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_qsc : PJ
	{
		FACE face;
		double a_squared;
		double b;
		double one_minus_f;
		double one_minus_f_squared;

		public override string Name { get { return "qsc"; } }
		public override string DescriptionName { get { return "Quadrilateralized Spherical Cube"; } }
		public override string DescriptionType { get { return "Azi, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double EPS10=1.0e-10;

		// The six cube faces.
		enum FACE
		{
			FRONT=0,
			RIGHT=1,
			BACK=2,
			LEFT=3,
			TOP=4,
			BOTTOM=5,
		}

		// The four areas on a cube face. AREA_0 is the area of definition,
		// the other three areas are counted counterclockwise.
		enum AREA
		{
			_0=0,
			_1=1,
			_2=2,
			_3=3,
		}

		// Helper function for forward projection: compute the theta angle and determine the area number.
		static double qsc_fwd_equat_face_theta(double phi, double y, double x, out AREA area)
		{
			if(phi<EPS10)
			{
				area=AREA._0;
				return 0.0;
			}

			double theta=Math.Atan2(y, x);
			if(Math.Abs(theta)<=Proj.FORTPI)
			{
				area=AREA._0;
			}
			else if(theta>Proj.FORTPI&&theta<=Proj.HALFPI+Proj.FORTPI)
			{
				area=AREA._1;
				theta-=Proj.HALFPI;
			}
			else if(theta>Proj.HALFPI+Proj.FORTPI||theta<=-(Proj.HALFPI+Proj.FORTPI))
			{
				area=AREA._2;
				theta=(theta>=0.0?theta-Proj.PI:theta+Proj.PI);
			}
			else
			{
				area=AREA._3;
				theta+=Proj.HALFPI;
			}

			return theta;
		}

		// Helper function: shift the longitude.
		static double qsc_shift_lon_origin(double lon, double offset)
		{
			double slon=lon+offset;
			if(slon<-Proj.PI) return slon+Proj.TWOPI;
			else if(slon>Proj.PI) return slon-Proj.TWOPI;
			return slon;
		}

		// Forward projection, ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			// Convert the geodetic latitude to a geocentric latitude.
			// This corresponds to the shift from the ellipsoid to the sphere
			// described in [LK12].
			double lat;
			if(es!=0) lat=Math.Atan(one_minus_f_squared*Math.Tan(lp.phi));
			else lat=lp.phi;

			// Convert the input lat, lon into theta, phi as used by QSC.
			// This depends on the cube face and the area on it.
			// For the top and bottom face, we can compute theta and phi
			// directly from phi, lam. For the other faces, we must use
			// unit sphere cartesian coordinates as an intermediate step.
			double lon=lp.lam;
			double theta, phi;
			AREA area;
			if(face==FACE.TOP)
			{
				phi=Proj.HALFPI-lat;
				if(lon>=Proj.FORTPI&&lon<=Proj.HALFPI+Proj.FORTPI)
				{
					area=AREA._0;
					theta=lon-Proj.HALFPI;
				}
				else if(lon>Proj.HALFPI+Proj.FORTPI||lon<=-(Proj.HALFPI+Proj.FORTPI))
				{
					area=AREA._1;
					theta=(lon>0.0?lon-Proj.PI:lon+Proj.PI);
				}
				else if(lon>-(Proj.HALFPI+Proj.FORTPI)&&lon<=-Proj.FORTPI)
				{
					area=AREA._2;
					theta=lon+Proj.HALFPI;
				}
				else
				{
					area=AREA._3;
					theta=lon;
				}
			}
			else if (face == FACE.BOTTOM)
			{
				phi = Proj.HALFPI + lat;
				if (lon >= Proj.FORTPI && lon <= Proj.HALFPI + Proj.FORTPI)
				{
					area = AREA._0;
					theta = -lon + Proj.HALFPI;
				}
				else if (lon < Proj.FORTPI && lon >= -Proj.FORTPI)
				{
					area = AREA._1;
					theta = -lon;
				}
				else if (lon < -Proj.FORTPI && lon >= -(Proj.HALFPI + Proj.FORTPI))
				{
					area = AREA._2;
					theta = -lon - Proj.HALFPI;
				}
				else
				{
					area = AREA._3;
					theta = (lon > 0.0 ? -lon + Proj.PI : -lon - Proj.PI);
				}
			}
			else
			{
				if (face == FACE.RIGHT) lon = qsc_shift_lon_origin(lon, +Proj.HALFPI);
				else if (face == FACE.BACK) lon = qsc_shift_lon_origin(lon, +Proj.PI);
				else if (face == FACE.LEFT) lon = qsc_shift_lon_origin(lon, -Proj.HALFPI);

				double sinlat = Math.Sin(lat);
				double coslat = Math.Cos(lat);
				double sinlon = Math.Sin(lon);
				double coslon = Math.Cos(lon);
				double q = coslat * coslon;
				double r = coslat * sinlon;
				double s = sinlat;

				if (face == FACE.FRONT)
				{
					phi = Math.Acos(q);
					theta = qsc_fwd_equat_face_theta(phi, s, r, out area);
				}
				else if (face == FACE.RIGHT)
				{
					phi = Math.Acos(r);
					theta = qsc_fwd_equat_face_theta(phi, s, -q, out area);
				}
				else if (face == FACE.BACK)
				{
					phi = Math.Acos(-q);
					theta = qsc_fwd_equat_face_theta(phi, s, -r, out area);
				}
				else if (face == FACE.LEFT)
				{
					phi = Math.Acos(-r);
					theta = qsc_fwd_equat_face_theta(phi, s, q, out area);
				}
				else
				{
					// Impossible
					phi = theta = 0.0;
					area = AREA._0;
				}
			}

			// Compute mu and nu for the area of definition.
			// For mu, see Eq. (3-21) in [OL76], but note the typos:
			// compare with Eq. (3-14). For nu, see Eq. (3-38).
			double mu=Math.Atan((12.0/Proj.PI)*(theta+Math.Acos(Math.Sin(theta)*Math.Cos(Proj.FORTPI))-Proj.HALFPI));
			double t=Math.Sqrt((1.0-Math.Cos(phi))/(Math.Cos(mu)*Math.Cos(mu))/(1.0-Math.Cos(Math.Atan(1.0/Math.Cos(theta)))));
			//double nu=atan(t); // We don't really need nu, just t, see below.

			// Apply the result to the real area.
			if(area==AREA._1) mu+=Proj.HALFPI;
			else if(area==AREA._2) mu+=Proj.PI;
			else if(area==AREA._3) mu+=Proj.HALFPI+Proj.PI;

			// Now compute x, y from mu and nu
			//t=tan(nu);
			xy.x=t*Math.Cos(mu);
			xy.y=t*Math.Sin(mu);

			return xy;
		}

		// Inverse projection, ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			// Convert the input x, y to the mu and nu angles as used by QSC.
			// This depends on the area of the cube face.
			double nu=Math.Atan(Math.Sqrt(xy.x*xy.x+xy.y*xy.y));
			double mu=Math.Atan2(xy.y, xy.x);
			AREA area;
			if(xy.x>=0.0&&xy.x>=Math.Abs(xy.y))
			{
				area=AREA._0;
			}
			else if(xy.y>=0.0&&xy.y>=Math.Abs(xy.x))
			{
				area=AREA._1;
				mu-=Proj.HALFPI;
			}
			else if(xy.x<0.0&&-xy.x>=Math.Abs(xy.y))
			{
				area=AREA._2;
				mu=(mu<0.0?mu+Proj.PI:mu-Proj.PI);
			}
			else
			{
				area=AREA._3;
				mu+=Proj.HALFPI;
			}

			// Compute phi and theta for the area of definition.
			// The inverse projection is not described in the original paper, but some
			// good hints can be found here (as of 2011-12-14):
			// http://fits.gsfc.nasa.gov/fitsbits/saf.93/saf.9302
			// (search for "Message-Id: <9302181759.AA25477 at fits.cv.nrao.edu>")
			double t=(Proj.PI/12.0)*Math.Tan(mu);
			double tantheta=Math.Sin(t)/(Math.Cos(t)-(1.0/Math.Sqrt(2.0)));
			double theta=Math.Atan(tantheta);
			double cosmu=Math.Cos(mu);
			double tannu=Math.Tan(nu);
			double cosphi=1.0-cosmu*cosmu*tannu*tannu*(1.0-Math.Cos(Math.Atan(1.0/Math.Cos(theta))));
			if(cosphi<-1.0) cosphi=-1.0;
			else if(cosphi>+1.0) cosphi=+1.0;

			// Apply the result to the real area on the cube face.
			// For the top and bottom face, we can compute phi and lam directly.
			// For the other faces, we must use unit sphere cartesian coordinates
			// as an intermediate step.
			double phi;

			if(face==FACE.TOP)
			{
				phi=Math.Acos(cosphi);
				lp.phi=Proj.HALFPI-phi;
				if(area==AREA._0) lp.lam=theta+Proj.HALFPI;
				else if(area==AREA._1) lp.lam=(theta<0.0?theta+Proj.PI:theta-Proj.PI);
				else if(area==AREA._2) lp.lam=theta-Proj.HALFPI;
				else lp.lam=theta; // area == AREA_3
			}
			else if(face==FACE.BOTTOM)
			{
				phi=Math.Acos(cosphi);
				lp.phi=phi-Proj.HALFPI;
				if(area==AREA._0) lp.lam=-theta+Proj.HALFPI;
				else if(area==AREA._1) lp.lam=-theta;
				else if(area==AREA._2) lp.lam=-theta-Proj.HALFPI;
				else lp.lam=(theta<0.0?-theta-Proj.PI:-theta+Proj.PI); // area == AREA_3
			}
			else
			{
				// Compute phi and lam via cartesian unit sphere coordinates.
				double q=cosphi;

				t=q*q;

				double s;
				if(t>=1.0) s=0.0;
				else s=Math.Sqrt(1.0-t)*Math.Sin(theta);

				t+=s*s;

				double r;
				if(t>=1.0) r=0.0;
				else r=Math.Sqrt(1.0-t);

				// Rotate q, r and s into the correct area.
				if(area==AREA._1)
				{
					t=r;
					r=-s;
					s=t;
				}
				else if(area==AREA._2)
				{
					r=-r;
					s=-s;
				}
				else if(area==AREA._3)
				{
					t=r;
					r=s;
					s=-t;
				}

				// Rotate q, r and s into the correct cube face.
				if(face==FACE.RIGHT)
				{
					t=q;
					q=-r;
					r=t;
				}
				else if(face==FACE.BACK)
				{
					q=-q;
					r=-r;
				}
				else if(face==FACE.LEFT)
				{
					t=q;
					q=r;
					r=-t;
				}

				// Now compute phi and lam from the unit sphere coordinates.
				lp.phi=Math.Acos(-s)-Proj.HALFPI;
				lp.lam=Math.Atan2(r, q);
				if(face==FACE.RIGHT) lp.lam=qsc_shift_lon_origin(lp.lam, -Proj.HALFPI);
				else if(face==FACE.BACK) lp.lam=qsc_shift_lon_origin(lp.lam, -Proj.PI);
				else if(face==FACE.LEFT) lp.lam=qsc_shift_lon_origin(lp.lam, +Proj.HALFPI);
			}

			// Apply the shift from the sphere to the ellipsoid as described in [LK12].
			if(es!=0)
			{
				bool invert_sign=lp.phi<0.0;

				double tanphi=Math.Tan(lp.phi);
				double xa=b/Math.Sqrt(tanphi*tanphi+one_minus_f_squared);
				lp.phi=Math.Atan(Math.Sqrt(a*a-xa*xa)/(one_minus_f*xa));
				if(invert_sign) lp.phi=-lp.phi;
			}

			return lp;
		}

		public override PJ Init()
		{
			fwd=e_forward;
			inv=e_inverse;

			// Determine the cube face from the center of projection.
			if(phi0>=Proj.HALFPI-Proj.FORTPI/2.0) face=FACE.TOP;
			else if(phi0<=-(Proj.HALFPI-Proj.FORTPI/2.0)) face=FACE.BOTTOM;
			else if(Math.Abs(lam0)<=Proj.FORTPI) face=FACE.FRONT;
			else if(Math.Abs(lam0)<=Proj.HALFPI+Proj.FORTPI) face=(lam0>0.0?FACE.RIGHT:FACE.LEFT);
			else face=FACE.BACK;

			// Fill in useful values for the ellipsoid <-> sphere shift
			// described in [LK12].
			if(es!=0)
			{
				a_squared=a*a;
				b=a*Math.Sqrt(1.0-es);
				one_minus_f=1.0-(a-b)/a;
				one_minus_f_squared=one_minus_f*one_minus_f;
			}

			return this;
		}
	}
}
