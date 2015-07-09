using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_putp4p : PJ
	{
		protected double C_x, C_y;

		public override string Name { get { return "putp4p"; } }
		public override string DescriptionName { get { return "Putnins P4'"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			lp.phi=Proj.aasin(ctx, 0.883883476*Math.Sin(lp.phi));
			xy.x=C_x*lp.lam*Math.Cos(lp.phi);
			lp.phi*=0.333333333333333;
			xy.x/=Math.Cos(lp.phi);
			xy.y=C_y*Math.Sin(lp.phi);

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			lp.phi=Proj.aasin(ctx, xy.y/C_y);
			lp.lam=xy.x*Math.Cos(lp.phi)/C_x;
			lp.phi*=3.0;
			lp.lam/=Math.Cos(lp.phi);
			lp.phi=Proj.aasin(ctx, 1.13137085*Math.Sin(lp.phi));

			return lp;
		}

		protected PJ_putp4p setup()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}

		public override PJ Init()
		{
			C_x=0.874038744;
			C_y=3.883251825;

			return setup();
		}
	}

	class PJ_weren : PJ_putp4p
	{
		public override string Name { get { return "weren"; } }
		public override string DescriptionName { get { return "Werenskiold I"; } }

		public override PJ Init()
		{
			C_x=1.0;
			C_y=4.442882938;

			return setup();
		}
	}
}
