// libproj -- library of cartographic projections
//
// Copyright (c) 2003, 2006 Gerald I. Evenden
// Copyright (c) 2009-2011 by the Authors
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

namespace Free.Ports.Proj4.Projections
{
	class PJ_rouss : PJ
	{
		protected double s0, A1, A2, A3, A4, A5, A6, B1, B2, B3, B4, B5, B6, B7, B8;
		protected double C1, C2, C3, C4, C5, C6, C7, C8, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11;
		protected MDIST en;

		public override string Name { get { return "rouss"; } }
		public override string DescriptionName { get { return "Roussilhe Stereographic"; } }
		public override string DescriptionType { get { return "Azi, Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double cp=Math.Cos(lp.phi);
			double sp=Math.Sin(lp.phi);
			double s=en.proj_mdist(lp.phi, sp, cp)-s0;
			double s2=s*s;
			double al=lp.lam*cp/Math.Sqrt(1.0-es*sp*sp);
			double al2=al*al;

			xy.x=k0*al*(1.0+s2*(A1+s2*A4)-al2*(A2+s*A3+s2*A5+al2*A6));
			xy.y=k0*(al2*(B1+al2*B4)+s*(1.0+al2*(B3-al2*B6)+s2*(B2+s2*B8)+s*al2*(B5+s*B7)));

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double x=xy.x/k0;
			double y=xy.y/k0;

			double x2=x*x;
			double y2=y*y;
			double al=x*(1.0-C1*y2+x2*(C2+C3*y-C4*x2+C5*y2-C7*x2*y)+y2*(C6*y2-C8*x2*y));
			double s=s0+y*(1.0+y2*(-D2+D8*y2))+x2*(-D1+y*(-D3+y*(-D5+y*(-D7+y*D11)))+x2*(D4+y*(D6+y*D10)-x2*D9));
			lp.phi=en.proj_inv_mdist(ctx, s);
			s=Math.Sin(lp.phi);
			lp.lam=al*Math.Sqrt(1.0-es*s*s)/Math.Cos(lp.phi);

			return lp;
		}

		public override PJ Init()
		{
			try { en=new MDIST(es); }
			catch { return null; }

			double es2=Math.Sin(phi0);
			s0=en.proj_mdist(phi0, es2, Math.Cos(phi0));
			es2=es*es2*es2;
			double t=1.0-es2;
			double N0=1.0/Math.Sqrt(t);
			double R_R0_2=t*t/one_es;
			double R_R0_4=R_R0_2*R_R0_2;
			t=Math.Tan(phi0);
			double t2=t*t;
			C1=A1=R_R0_2/4.0;
			C2=A2=R_R0_2*(2*t2-1.0-2.0*es2)/12.0;
			A3=R_R0_2*t*(1.0+4.0*t2)/(12.0*N0);
			A4=R_R0_4/24.0;
			A5=R_R0_4*(-1.0+t2*(11.0+12.0*t2))/24.0;
			A6=R_R0_4*(-2.0+t2*(11.0-2.0*t2))/240.0;
			B1=t/(2.0*N0);
			B2=R_R0_2/12.0;
			B3=R_R0_2*(1.0+2.0*t2-2.0*es2)/4.0;
			B4=R_R0_2*t*(2.0-t2)/(24.0*N0);
			B5=R_R0_2*t*(5.0+4.0*t2)/(8.0*N0);
			B6=R_R0_4*(-2.0+t2*(-5.0+6.0*t2))/48.0;
			B7=R_R0_4*(5.0+t2*(19.0+12.0*t2))/24.0;
			B8=R_R0_4/120.0;
			C3=R_R0_2*t*(1.0+t2)/(3.0*N0);
			C4=R_R0_4*(-3.0+t2*(34.0+22.0*t2))/240.0;
			C5=R_R0_4*(4.0+t2*(13.0+12.0*t2))/24.0;
			C6=R_R0_4/16.0;
			C7=R_R0_4*t*(11.0+t2*(33.0+t2*16.0))/(48.0*N0);
			C8=R_R0_4*t*(1.0+t2*4.0)/(36.0*N0);
			D1=t/(2.0*N0);
			D2=R_R0_2/12.0;
			D3=R_R0_2*(2*t2+1.0-2.0*es2)/4.0;
			D4=R_R0_2*t*(1.0+t2)/(8.0*N0);
			D5=R_R0_2*t*(1.0+t2*2.0)/(4.0*N0);
			D6=R_R0_4*(1.0+t2*(6.0+t2*6.0))/16.0;
			D7=R_R0_4*t2*(3.0+t2*4.0)/8.0;
			D8=R_R0_4/80.0;
			D9=R_R0_4*t*(-21.0+t2*(178.0-t2*26.0))/720.0;
			D10=R_R0_4*t*(29.0+t2*(86.0+t2*48.0))/(96.0*N0);
			D11=R_R0_4*t*(37.0+t2*44.0)/(96.0*N0);
			fwd=e_forward;
			inv=e_inverse;

			return this;
		}
	}
}
