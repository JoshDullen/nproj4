using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_putp6 : PJ
	{
		protected double C_x, C_y, A, B, D;

		public override string Name { get { return "putp6"; } }
		public override string DescriptionName { get { return "Putnins P6"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double EPS=1.0e-10;
		const int NITER=10;
		const double CON_POLE=1.732050807568877;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double p=B*Math.Sin(lp.phi);
			lp.phi*=1.10265779;

			int i=NITER;
			for(; i>0; i--)
			{
				double r=Math.Sqrt(1.0+lp.phi*lp.phi);
				double V=((A-r)*lp.phi-Math.Log(lp.phi+r)-p)/(A-2.0*r);
				lp.phi-=V;
				if(Math.Abs(V)<EPS) break;
			}

			if(i==0) lp.phi=p<0.0?-CON_POLE:CON_POLE;

			xy.x=C_x*lp.lam*(D-Math.Sqrt(1.0+lp.phi*lp.phi));
			xy.y=C_y*lp.phi;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y/C_y;
			double r=Math.Sqrt(1.0+lp.phi*lp.phi);
			lp.lam=xy.x/(C_x*(D-r));
			lp.phi=Proj.aasin(ctx, ((A-r)*lp.phi-Math.Log(lp.phi+r))/B);

			return lp;
		}

		protected PJ_putp6 setup()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}

		public override PJ Init()
		{
			C_x=1.01346;
			C_y=0.91910;
			A=4.0;
			B=2.1471437182129378784;
			D=2.0;

			return setup();
		}
	}

	class PJ_putp6p : PJ_putp6
	{
		public override string Name { get { return "putp6p"; } }
		public override string DescriptionName { get { return "Putnins P6'"; } }

		public override PJ Init()
		{
			C_x=0.44329;
			C_y=0.80404;
			A=6.0;
			B=5.61125;
			D=3.0;

			return setup();
		}
	}
}
