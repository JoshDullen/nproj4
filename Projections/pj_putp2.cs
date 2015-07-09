using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_putp2 : PJ
	{
		public override string Name { get { return "putp2"; } }
		public override string DescriptionName { get { return "Putnins P2"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double C_x=1.89490;
		const double C_y=1.71848;
		const double C_p=0.6141848493043784;
		const double EPS=1.0e-10;
		const int NITER=10;
		const double PI_DIV_3=1.0471975511965977;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double p=C_p*Math.Sin(lp.phi);
			double s=lp.phi*lp.phi;
			lp.phi*=0.615709+s*(0.00909953+s*0.0046292);

			int i=NITER;
			for(; i>0; i--)
			{
				double c=Math.Cos(lp.phi);
				s=Math.Sin(lp.phi);
				double V=(lp.phi+s*(c-1.0)-p)/(1.0+c*(c-1.0)-s*s);
				lp.phi-=V;
				if(Math.Abs(V)<EPS) break;
			}

			if(i==0) lp.phi=lp.phi<0?-PI_DIV_3:PI_DIV_3;

			xy.x=C_x*lp.lam*(Math.Cos(lp.phi)-0.5);
			xy.y=C_y*Math.Sin(lp.phi);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=Proj.aasin(ctx, xy.y/C_y);
			double c=Math.Cos(lp.phi);
			lp.lam=xy.x/(C_x*(c-0.5));
			lp.phi=Proj.aasin(ctx, (lp.phi+Math.Sin(lp.phi)*(c-1.0))/C_p);

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
