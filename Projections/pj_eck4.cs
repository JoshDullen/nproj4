using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_eck4 : PJ
	{
		public override string Name { get { return "eck4"; } }
		public override string DescriptionName { get { return "Eckert IV"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double C_x=0.42223820031577120149;
		const double C_y=1.32650042817700232218;
		const double RC_y=0.75386330736002178205;
		const double C_p=3.57079632679489661922;
		const double RC_p=0.28004957675577868795;
		const double EPS=1.0e-7;
		const int NITER=6;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double p=C_p*Math.Sin(lp.phi);
			double V=lp.phi*lp.phi;
			lp.phi*=0.895168+V*(0.0218849+V*0.00826809);

			int i=NITER;
			for(; i>0; i--)
			{
				double c=Math.Cos(lp.phi);
				double s=Math.Sin(lp.phi);
				lp.phi-=V=(lp.phi+s*(c+2.0)-p)/(1.0+c*(c+2.0)-s*s);
				if(Math.Abs(V)<EPS) break;
			}

			if(i==0)
			{
				xy.x=C_x*lp.lam;
				xy.y=lp.phi<0.0?-C_y:C_y;
			}
			else
			{
				xy.x=C_x*lp.lam*(1.0+Math.Cos(lp.phi));
				xy.y=C_y*Math.Sin(lp.phi);
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=Proj.aasin(ctx, xy.y/C_y);
			double c=Math.Cos(lp.phi);
			lp.lam=xy.x/(C_x*(1.0+c));
			lp.phi=Proj.aasin(ctx, (lp.phi+Math.Sin(lp.phi)*(c+2.0))/C_p);

			return lp;
		}

		public override PJ Init()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}
