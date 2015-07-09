using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	// libproj -- library of cartographic projections
	//
	// Copyright (c) 2006-2008 Gerald I. Evenden
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
	class PJ_ftmerc : PJ
	{
		protected double rho0;
		protected double[] fc=new double[5], ic=new double[5];

		public override string Name { get { return "ftmerc"; } }
		public override string DescriptionName { get { return "French Transverse Mercator"; } }
		public override string DescriptionType { get { return "Cyl, Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// forward series constants
		const double FC00=1.0;
		const double FC02=-0.25;
		const double FC04=-0.046875;
		const double FC06=-0.01953125;
		const double FC08=-0.01068115234375;
		const double FC22=0.125;
		const double FC24=-0.01041666666666666666666666667;
		const double FC26=-0.0087890625;
		const double FC28=-0.004888237847222222222222222222;
		const double FC44=0.01692708333333333333333333333;
		const double FC46=0.0033203125;
		const double FC48=-0.0004218207465277777777777777778;
		const double FC66=0.003971354166666666666666666667;
		const double FC68=0.002090308779761904761904761905;
		const double FC88=0.001200382293216765873015873016;

		// inverse series constants
		const double IC00=1.0;
		const double IC02=-0.25;
		const double IC04=-0.046875;
		const double IC06=-0.01953125;
		const double IC08=-0.01068115234375;
		const double IC22=0.125;
		const double IC24=0.02083333333333333333333333333;
		const double IC26=0.00341796875;
		const double IC28=0.00001627604166666666666666666667;
		const double IC44=0.001302083333333333333333333333;
		const double IC46=0.00234375;
		const double IC48=0.001516384548611111111111111111;
		const double IC66=0.0005533854166666666666666666667;
		const double IC68=0.0006580171130952380952380952381;
		const double IC88=0.0001064966595362103174603174603;

		static void fs_init(double[] s, double es)
		{
			double t=es*es;
			s[0]=FC00+es*(FC02+es*(FC04+es*(FC06+es*FC08)));
			s[1]=es*(FC22+es*(FC24+es*(FC26+es*FC28)));
			s[2]=t*(FC44+es*(FC46+es*FC48));
			t*=es;
			s[3]=t*(FC66+es*FC68);
			s[4]=t*es*FC88;
		}

		static void is_init(double[] s, double es)
		{
			double t=es*es;
			s[0]=IC00+es*(IC02+es*(IC04+es*(IC06+es*IC08)));
			s[1]=es*(IC22+es*(IC24+es*(IC26+es*IC28)));
			s[2]=t*(IC44+es*(IC46+es*IC48));
			t*=es;
			s[3]=t*(IC66+es*IC68);
			s[4]=t*es*IC88;
		}

		static Complex cevals(Complex z, double[] c)
		{
			int ci=0;
			Complex zp=z*c[ci];
			for(int i=2; i<=8; i+=2) zp+=c[++ci]*Complex.Sin(i*z);
			return zp;
		}

		static Complex icevals(Complex z, double[] c)
		{
			int ci=0;
			Complex zp=z/=c[ci];
			for(int i=2; i<=8; i+=2) zp-=c[++ci]*Complex.Sin(i*z);
			return zp;
		}

		// isometric latitude
		static double proj_psi(double phi, double sphi, double e)
		{
			double esp=e*sphi;

			return Math.Log(Math.Tan(Proj.FORTPI+0.5*phi)*Math.Pow((1.0-esp)/(1.0+esp), 0.5*e));
		}

		// inverse isometric latitude
		static double proj_apsi(double psi, double e)
		{
			const int MAX_ITER=11;
			const double EPS=1e-14;

			double he=e*0.5, exp_psi=Math.Exp(psi);
			int i=MAX_ITER;

			double phi=0.0;
			double phi0=2.0*Math.Atan(exp_psi)-Proj.HALFPI;
			while((Math.Abs(phi0-phi)>EPS)&&(--i)>0)
			{
				double esp=e*Math.Sin(phi0);
				phi=2.0*Math.Atan(Math.Pow((1.0+esp)/(1.0-esp), he)*exp_psi)-Proj.HALFPI;
				phi0=phi;
			}

			// if(!i) runaway
			return phi;
		}

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double psi=proj_psi(lp.phi, Math.Sin(lp.phi), e);
			double beta=Math.Asin(Math.Sin(lp.lam)/Math.Cosh(psi));
			double psi_s=Math.Log(Math.Tan(Proj.FORTPI+0.5*beta));
			Complex z=new Complex(Math.Atan(Math.Sinh(psi)/Math.Cos(lp.lam)), psi_s);
			z=cevals(z, fc);
			xy.x=z.im*k0;
			xy.y=(z.re-rho0)*k0;

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.x/=k0;
			xy.y=xy.y/k0+rho0;
			Complex z=new Complex(xy.y, xy.x);
			z=icevals(z, ic);
			double L=z.re;
			double Ls=z.im;
			lp.lam=Math.Atan(Math.Sinh(Ls)/Math.Cos(L));
			double psi=Math.Asin(Math.Sin(L)/Math.Cosh(Ls));
			L=Math.Log(Math.Tan(Proj.FORTPI+0.5*psi));
			lp.phi=proj_apsi(L, e);

			return lp;
		}

		public override PJ Init()
		{
			fs_init(fc, es);
			is_init(ic, es);
			try { rho0=new MDIST(es).proj_mdist(phi0, Math.Sin(phi0), Math.Cos(phi0)); }
			catch { }
			fwd=e_forward;
			inv=e_inverse;

			return this;
		}

		// Revision 3.4 2008/06/26 15:18:06 gie
		// some minor repairs, still non-C90
		//
		// Revision 3.3 2006/06/19 01:06:04 gie
		// removed 'dummy' from entry
		//
		// Revision 3.2 2006/06/19 00:58:58 gie
		// fix Id by adding colon
		//
		// Revision 3.1 2006/06/19 00:54:57 gie
		// initial
	}
}
