using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_mbt_fps : PJ
	{
		public override string Name { get { return "mbt_fps"; } }
		public override string DescriptionName { get { return "McBryde-Thomas Flat-Pole Sine (No. 2)"; } }
		public override string DescriptionType { get { return "Cyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const int MAX_ITER=10;
		const double LOOP_TOL=1.0e-7;

		const double C1=0.45503;
		const double C2=1.36509;
		const double C3=1.41546;
		const double C_x=0.22248;
		const double C_y=1.44492;
		const double C1_2=0.33333333333333333333333333;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double k=C3*Math.Sin(lp.phi);
			double t;

			for(int i=MAX_ITER; i>0; i--)
			{
				t=lp.phi/C2;
				double V=(C1*Math.Sin(t)+Math.Sin(lp.phi)-k)/(C1_2*Math.Cos(t)+Math.Cos(lp.phi));
				lp.phi-=V;
				if(Math.Abs(V)<LOOP_TOL) break;
			}

			t=lp.phi/C2;
			xy.x=C_x*lp.lam*(1.0+3.0*Math.Cos(lp.phi)/Math.Cos(t));
			xy.y=C_y*Math.Sin(t);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double t=Proj.aasin(ctx, xy.y/C_y);

			lp.phi=C2*t;
			lp.lam=xy.x/(C_x*(1.0+3.0*Math.Cos(lp.phi)/Math.Cos(t)));
			lp.phi=Proj.aasin(ctx, (C1*Math.Sin(t)+Math.Sin(lp.phi))/C3);

			return lp;
		}

		public override PJ Init()
		{
			es=0;
			inv=s_inverse;
			fwd=s_forward;
			
			return this;
		}
	}
}
