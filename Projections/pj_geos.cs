//
// libproj -- library of cartographic projections
//
// Copyright (c) 2004 Gerald I. Evenden
// Copyright (c) 2012 Martin Raspaud
// Copyright (c) 2008-2015 by the Authors
//
// See also (section 4.4.3.2):
//	http://www.eumetsat.int/en/area4/msg/news/us_doc/cgms_03_26.pdf
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
// CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
// SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

using System;
using Free.Ports.Proj4.LibcStuff;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_geos : PJ
	{
		protected double h, radius_p, radius_p2, radius_p_inv2, radius_g, radius_g_1, C;
		protected bool flip_axis;

		public override string Name { get { return "geos"; } }
		public override string DescriptionName { get { return "Geostationary Satellite View"; } }
		public override string DescriptionType { get { return "Azi, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "h= sweep="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +h={0}", h);
				if(flip_axis) ret.Append(" +sweep=x");
				return ret.ToString();
			}
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			// Calculation of the three components of the vector from satellite to
			// position on earth surface (lon,lat).
			double tmp=Math.Cos(lp.phi);
			double Vx=Math.Cos(lp.lam)*tmp;
			double Vy=Math.Sin(lp.lam)*tmp;
			double Vz=Math.Sin(lp.phi);

			// Check visibility.
			if(((radius_g-Vx)*Vx-Vy*Vy-Vz*Vz)<0.0) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			// Calculation based on view angles from satellite.
			tmp=radius_g-Vx;
			if(flip_axis)
			{
				xy.x=radius_g_1*Math.Atan(Vy/Libc.hypot(Vz, tmp));
				xy.y=radius_g_1*Math.Atan(Vz/tmp);
			}
			else
			{
				xy.x=radius_g_1*Math.Atan(Vy/tmp);
				xy.y=radius_g_1*Math.Atan(Vz/Libc.hypot(Vy, tmp));
			}

			return xy;
		}

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			// Calculation of geocentric latitude.
			lp.phi=Math.Atan(radius_p2*Math.Tan(lp.phi));

			// Calculation of the three components of the vector from satellite to
			// position on earth surface (lon,lat).
			double r=radius_p/Libc.hypot(radius_p*Math.Cos(lp.phi), Math.Sin(lp.phi));
			double Vx=r*Math.Cos(lp.lam)*Math.Cos(lp.phi);
			double Vy=r*Math.Sin(lp.lam)*Math.Cos(lp.phi);
			double Vz=r*Math.Sin(lp.phi);

			// Check visibility.
			if(((radius_g-Vx)*Vx-Vy*Vy-Vz*Vz*radius_p_inv2)<0.0) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			// Calculation based on view angles from satellite.
			double tmp=radius_g-Vx;
			if(flip_axis)
			{
				xy.x=radius_g_1*Math.Atan(Vy/Libc.hypot(Vz, tmp));
				xy.y=radius_g_1*Math.Atan(Vz/tmp);
			}
			else
			{
				xy.x=radius_g_1*Math.Atan(Vy/tmp);
				xy.y=radius_g_1*Math.Atan(Vz/Libc.hypot(Vy, tmp));
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			// Setting three components of vector from satellite to position.
			double Vx=-1.0;
			double Vy, Vz;
			if(flip_axis)
			{
				Vz=Math.Tan(xy.y/(radius_g-1.0));
				Vy=Math.Tan(xy.x/(radius_g-1.0))*Math.Sqrt(1.0+Vz*Vz);
			}
			else
			{
				Vy=Math.Tan(xy.x/(radius_g-1.0));
				Vz=Math.Tan(xy.y/(radius_g-1.0))*Math.Sqrt(1.0+Vy*Vy);
			}

			// Calculation of terms in cubic equation and determinant.
			double a=Vy*Vy+Vz*Vz+Vx*Vx;
			double b=2*radius_g*Vx;
			double det=(b*b)-4*a*C;
			if(det<0.0) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

			// Calculation of three components of vector from satellite to position.
			double k=(-b-Math.Sqrt(det))/(2*a);
			Vx=radius_g+k*Vx;
			Vy*=k;
			Vz*=k;

			// Calculation of longitude and latitude.
			lp.lam=Math.Atan2(Vy, Vx);
			lp.phi=Math.Atan(Vz*Math.Cos(lp.lam)/Vx);

			return lp;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			// Setting three components of vector from satellite to position.
			double Vx=-1.0;
			double Vy, Vz;
			if(flip_axis)
			{
				Vz=Math.Tan(xy.y/radius_g_1);
				Vy=Math.Tan(xy.x/radius_g_1)*Libc.hypot(1.0, Vz);
			}
			else
			{
				Vy=Math.Tan(xy.x/radius_g_1);
				Vz=Math.Tan(xy.y/radius_g_1)*Libc.hypot(1.0, Vy);
			}

			// Calculation of terms in cubic equation and determinant.
			double a=Vz/radius_p;
			a=Vy*Vy+a*a+Vx*Vx;
			double b=2*radius_g*Vx;
			double det=(b*b)-4*a*C;
			if(det<0.0) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

			// Calculation of three components of vector from satellite to position.
			double k=(-b-Math.Sqrt(det))/(2.0*a);
			Vx=radius_g+k*Vx;
			Vy*=k;
			Vz*=k;

			// Calculation of longitude and latitude.
			lp.lam=Math.Atan2(Vy, Vx);
			lp.phi=Math.Atan(Vz*Math.Cos(lp.lam)/Vx);
			lp.phi=Math.Atan(radius_p_inv2*Math.Tan(lp.phi));

			return lp;
		}

		public override PJ Init()
		{
			h=Proj.pj_param_d(ctx, parameters, "h");
			if(h<=0.0) { Proj.pj_ctx_set_errno(ctx, -30); return null; }
			if(phi0!=0) { Proj.pj_ctx_set_errno(ctx, -46); return null; }
			string sweep_axis=Proj.pj_param_s(ctx, parameters, "sweep");
			if(sweep_axis==null||sweep_axis=="") flip_axis=false;
			else
			{
				if(sweep_axis.Length>1||(sweep_axis[0]!='x'&&sweep_axis[0]!='y'))
				{ Proj.pj_ctx_set_errno(ctx, -49); return null; }

				flip_axis=sweep_axis[0]=='x';
			}
			radius_g_1=h/a;
			radius_g=1.0+radius_g_1;
			C=radius_g*radius_g-1.0;

			if(es!=0)
			{
				radius_p=Math.Sqrt(one_es);
				radius_p2=one_es;
				radius_p_inv2=rone_es;
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{
				radius_p=radius_p2=radius_p_inv2=1.0;
				inv=s_inverse;
				fwd=s_forward;
			}

			return this;
		}
	}
}
