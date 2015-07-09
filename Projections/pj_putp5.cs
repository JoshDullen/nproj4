using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_putp5 : PJ
	{
		protected double A, B;

		public override string Name { get { return "putp5"; } }
		public override string DescriptionName { get { return "Putnins P5"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double C=1.01346;
		const double D=1.2158542;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=C*lp.lam*(A-B*Math.Sqrt(1.0+D*lp.phi*lp.phi));
			xy.y=C*lp.phi;

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=xy.y/C;
			lp.lam=xy.x/(C*(A-B*Math.Sqrt(1.0+D*lp.phi*lp.phi)));

			return lp;
		}

		protected PJ_putp5 setup()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}

		public override PJ Init()
		{
			A=2.0;
			B=1.0;

			return setup();
		}
	}

	class PJ_putp5p : PJ_putp5
	{
		public override string Name { get { return "putp5p"; } }
		public override string DescriptionName { get { return "Putnins P5'"; } }

		public override PJ Init()
		{
			A=1.5;
			B=0.5;

			return setup();
		}
	}
}
