using System;

namespace Free.Ports.Proj4.Projections
{
	// libproj -- library of cartographic projections
	//
	// Copyright (c) 2008   Gerald I. Evenden
	// Copyright (c) 2008-2011 by the Authors
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
	class PJ_ktmerc : PJ
	{
		protected double As, Bs, Cs, Ds; // conformal lat constants
		protected double Ap, Bp, Cp, Dp; // conformal lat constants
		protected double[] beta=new double[4], del=new double[4];
		protected double K;

		public override string Name { get { return "ktmerc"; } }
		public override string DescriptionName { get { return "Kruger Transverse Mercator"; } }
		public override string DescriptionType { get { return "Cyl, Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		public static double Atanh(double z)
		{
			return 0.5*Math.Log((1+z)/(1-z));
		}

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sphi2=Math.Sin(lp.phi);
			double c=Math.Cos(lp.phi)*sphi2;
			sphi2*=sphi2;
			double phis=lp.phi-c*(As+sphi2*(Bs+sphi2*(Cs+sphi2*Ds)));
			double z=xy.y=Math.Atan2(Math.Tan(phis), Math.Cos(lp.lam));
			double n=xy.x=Atanh(Math.Cos(phis)*Math.Sin(lp.lam));
			for(int b=0, i=2; i<=8; i+=2, ++b)
			{
				double tz=i*z;
				double tn=i*n;
				xy.y+=beta[b]*Math.Sin(tz)*Math.Cosh(tn);
				xy.x+=beta[b]*Math.Cos(tz)*Math.Sinh(tn);
			}
			xy.y*=K;
			xy.x*=K;
			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double z, n;
			double sz=z=xy.y/K;
			double sn=n=xy.x/K;
			for(int d=0, i=2; i<=8; i+=2, ++d)
			{
				double tz=i*z;
				double tn=i*n;
				sz-=del[d]*Math.Sin(tz)*Math.Cosh(tn);
				sn-=del[d]*Math.Cos(tz)*Math.Sinh(tn);
			}

			lp.lam=Math.Atan2(Math.Sinh(sn), Math.Cos(sz));
			lp.phi=Math.Asin(Math.Sin(sz)/Math.Cosh(sn));
			double sinp2=Math.Sin(lp.phi);
			double c=Math.Cos(lp.phi)*sinp2;
			sinp2*=sinp2;
			lp.phi+=c*(Ap+sinp2*(Bp+sinp2*(Cp+sinp2*Dp)));

			return lp;
		}

		public override PJ Init()
		{
			if(es<=0.0) { Proj.pj_ctx_set_errno(ctx, -34); return null; }
			double f=1.0-Math.Sqrt(one_es);
			double es4=es*es;
			double es6=es4*es;
			double es8=es6*es;
			As=es;
			Bs=es4*(5.0-es)/6.0;
			Cs=(es6*104.0-es8*45.0)/120.0;
			Ds=es8*1237.0/1260.0;
			Ap=es+es4+es6+es8;
			Bp=-(es4*7.0+es6*17.0+es8*30.0)/6.0;
			Cp=(es6*224.0+es8*889.0)/120;
			Dp=4279.0*es8/1260.0;

			double n=f/(2.0-f);
			double n2=n*n;
			K=k0*(1.0+n2*(1.0/4.0+n2*1.0/64.0))/(1+n);
			beta[0]=n*(1.0/2.0+n*(-2.0/3.0+n*(5.0/16.0+n*41.0/180.0)));
			beta[1]=n2*(13.0/48.0+n*(-3.0/5.0+n*557.0/1440.0));
			beta[2]=n2*n*(61.0/240.0-n*103.0/140.0);
			beta[3]=n2*n2*49561.0/161280.0;
			del[0]=n*(1.0/2.0+n*(-2.0/3.0+n*(37.0/96.0-n*1.0/360.0)));
			del[1]=n2*(1.0/48.0+n*(1.0/15.0-n*437.0/1440.0));
			del[2]=n*n2*(17.0/480.0-37.0/840.0);
			del[3]=n2*n2*4397.0/161280.0;

			inv=e_inverse;
			fwd=e_forward;

			return this;
		}
	}
}
