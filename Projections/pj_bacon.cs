using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_bacon : PJ
	{
		protected bool bacn;
		protected bool ortl;

		public override string Name { get { return "bacon"; } }
		public override string DescriptionName { get { return "Bacon Globular"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return false; } }

		const double HLFPI2=2.46740110027233965467;
		const double EPS=1.0e-10;

		// spheroid
		protected XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.y=bacn?Proj.HALFPI*Math.Sin(lp.phi):lp.phi;
			double ax=Math.Abs(lp.lam);

			if(ax>=EPS)
			{
				if(ortl&&ax>=Proj.HALFPI) xy.x=Math.Sqrt(HLFPI2-lp.phi*lp.phi+EPS)+ax-Proj.HALFPI;
				else
				{
					double f=0.5*(HLFPI2/ax+ax);
					xy.x=ax-f+Math.Sqrt(f*f-xy.y*xy.y);
				}
				if(lp.lam<0.0) xy.x=-xy.x;
			}
			else xy.x=0.0;

			return xy;
		}

		public override PJ Init()
		{
			bacn=true;
			ortl=false;
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}

	class PJ_apian : PJ_bacon
	{
		public override string Name { get { return "apian"; } }
		public override string DescriptionName { get { return "Apian Globular I"; } }

		public override PJ Init()
		{
			bacn=ortl=false;
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}

	class PJ_ortel : PJ_bacon
	{
		public override string Name { get { return "ortel"; } }
		public override string DescriptionName { get { return "Ortelius Oval"; } }

		public override PJ Init()
		{
			bacn=false;
			ortl=true;
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}
}
