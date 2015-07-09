using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_putp3 : PJ
	{
		protected double A;

		public override string Name { get { return "putp3"; } }
		public override string DescriptionName { get { return "Putnins P3"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double C=0.79788456;
		protected const double RPISQ=0.1013211836;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=C*lp.lam*(1.0-A*lp.phi*lp.phi);
			xy.y=C*lp.phi;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y/C;
			lp.lam=xy.x/(C*(1.0-A*lp.phi*lp.phi));

			return lp;
		}

		protected PJ_putp3 setup()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}

		public override PJ Init()
		{
			A=4.0*RPISQ;

			return setup();
		}
	}

	class PJ_putp3p : PJ_putp3
	{
		public override string Name { get { return "putp3p"; } }
		public override string DescriptionName { get { return "Putnins P3'"; } }

		public override PJ Init()
		{
			A=2.0*RPISQ;

			return setup();
		}
	}
}
