using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_moll : PJ
	{
		protected double C_x, C_y, C_p;

		public override string Name { get { return "moll"; } }
		public override string DescriptionName { get { return "Mollweide"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const int MAX_ITER=10;
		const double LOOP_TOL=1.0e-7;

		// spheroid
		protected XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double k=C_p*Math.Sin(lp.phi);
			int i=MAX_ITER;
			for(; i>0; i--)
			{
				double V=(lp.phi+Math.Sin(lp.phi)-k)/(1.0+Math.Cos(lp.phi));
				lp.phi-=V;
				if(Math.Abs(V)<LOOP_TOL) break;
			}

			if(i==0) lp.phi=(lp.phi<0.0)?-Proj.HALFPI:Proj.HALFPI;
			else lp.phi*=0.5;

			xy.x=C_x*lp.lam*Math.Cos(lp.phi);
			xy.y=C_y*Math.Sin(lp.phi);

			return xy;
		}

		// spheroid
		protected LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=Proj.aasin(ctx, xy.y/C_y);
			lp.lam=xy.x/(C_x*Math.Cos(lp.phi));
			lp.phi+=lp.phi;
			lp.phi=Proj.aasin(ctx, (lp.phi+Math.Sin(lp.phi))/C_p);

			return lp;
		}

		protected PJ_moll setup(double p)
		{
			double p2=p+p;
			es=0;
			double sp=Math.Sin(p);
			double r=Math.Sqrt(Proj.TWOPI*sp/(p2+Math.Sin(p2)));
			C_x=2.0*r/Proj.PI;
			C_y=r/sp;
			C_p=p2+Math.Sin(p2);
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}

		public override PJ Init()
		{
			return setup(Proj.HALFPI);
		}
	}

	class PJ_wag4 : PJ_moll
	{
		public override string Name { get { return "wag4"; } }
		public override string DescriptionName { get { return "Wagner IV"; } }

		public override PJ Init()
		{
			return setup(Proj.PI/3.0);
		}
	}

	class PJ_wag5 : PJ_moll
	{
		public override string Name { get { return "wag5"; } }
		public override string DescriptionName { get { return "Wagner V"; } }

		public override PJ Init()
		{
			es=0;
			C_x=0.90977;
			C_y=1.65014;
			C_p=3.00896;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}
