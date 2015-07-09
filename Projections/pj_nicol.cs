using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_nicol : PJ
	{
		public override string Name { get { return "nicol"; } }
		public override string DescriptionName { get { return "Nicolosi Globular"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const double EPS=1.0e-10;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(Math.Abs(lp.lam)<EPS)
			{
				xy.x=0;
				xy.y=lp.phi;
			}
			else if(Math.Abs(lp.phi)<EPS)
			{
				xy.x=lp.lam;
				xy.y=0.0;
			}
			else if(Math.Abs(Math.Abs(lp.lam)-Proj.HALFPI)<EPS)
			{
				xy.x=lp.lam*Math.Cos(lp.phi);
				xy.y=Proj.HALFPI*Math.Sin(lp.phi);
			}
			else if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<EPS)
			{
				xy.x=0;
				xy.y=lp.phi;
			}
			else
			{
				double tb=Proj.HALFPI/lp.lam-lp.lam/Proj.HALFPI;
				double c=lp.phi/Proj.HALFPI;
				double sp=Math.Sin(lp.phi);
				double d=(1-c*c)/(sp-c);
				double r2=tb/d;
				r2*=r2;
				double m=(tb*sp/d-0.5*tb)/(1.0+r2);
				double n=(sp/r2+0.5*d)/(1.0+1.0/r2);
				xy.x=Math.Cos(lp.phi);
				xy.x=Math.Sqrt(m*m+xy.x*xy.x/(1.0+r2));
				xy.x=Proj.HALFPI*(m+(lp.lam<0.0?-xy.x:xy.x));
				xy.y=Math.Sqrt(n*n-(sp*sp/r2+d*sp-1.0)/(1.0+1.0/r2));
				xy.y=Proj.HALFPI*(n+(lp.phi<0.0?xy.y:-xy.y));
			}

			return xy;
		}

		public override PJ Init()
		{
			es=0;
			fwd=s_forward;

			return this;
		}
	}
}
