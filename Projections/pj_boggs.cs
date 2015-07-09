using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_boggs : PJ
	{
		public override string Name { get { return "boggs"; } }
		public override string DescriptionName { get { return "Boggs Eumorphic"; } }
		public override string DescriptionType { get { return "PCyl, Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const int NITER=20;
		const double EPS=1.0e-7;
		const double FXC=2.00276;
		const double FXC2=1.11072;
		const double FYC=0.49931;
		const double FYC2=1.41421356237309504880;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double theta=lp.phi;
			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)>=EPS)
			{
				double c=Math.Sin(theta)*Proj.PI;
				for(int i=NITER; i>0; i--)
				{
					double th1=(theta+Math.Sin(theta)-c)/(1.0+Math.Cos(theta));
					theta-=th1;
					if(Math.Abs(th1)<EPS) break;
				}
				theta*=0.5;
				xy.x=FXC*lp.lam/(1.0/Math.Cos(lp.phi)+FXC2/Math.Cos(theta));
			}

			xy.y=FYC*(lp.phi+FYC2*Math.Sin(theta));

			return xy;
		}

		public override PJ Init()
		{
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}
}
