using System;

namespace Free.Ports.Proj4.Projections
{
	abstract class PJ_sts : PJ
	{
		protected double C_x, C_y, C_p;
		protected bool tan_mode;

		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			xy.x=C_x*lp.lam*Math.Cos(lp.phi);
			xy.y=C_y;
			lp.phi*=C_p;
			double c=Math.Cos(lp.phi);

			if(tan_mode)
			{
				xy.x*=c*c;
				xy.y*=Math.Tan(lp.phi);
			}
			else
			{
				xy.x/=c;
				xy.y*=Math.Sin(lp.phi);
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y/=C_y;
			lp.phi=tan_mode?Math.Atan(xy.y):Proj.aasin(ctx, xy.y);
			double c=Math.Cos(lp.phi);
			lp.phi/=C_p;
			lp.lam=xy.x/(C_x*Math.Cos(lp.phi));
			if(tan_mode) lp.lam/=c*c;
			else lp.lam*=c;

			return lp;
		}

		protected PJ_sts setup(double p, double q, bool mode)
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;
			C_x=q/p;
			C_y=p;
			C_p=1/q;
			tan_mode=mode;

			return this;
		}
	}

	class PJ_kav5 : PJ_sts
	{
		public override string Name { get { return "kav5"; } }
		public override string DescriptionName { get { return "Kavraisky V"; } }

		public override PJ Init()
		{
			return setup(1.50488, 1.35439, false);
		}
	}

	class PJ_qua_aut : PJ_sts
	{
		public override string Name { get { return "qua_aut"; } }
		public override string DescriptionName { get { return "Quartic Authalic"; } }

		public override PJ Init()
		{
			return setup(2, 2, false);
		}
	}

	class PJ_mbt_s : PJ_sts
	{
		public override string Name { get { return "mbt_s"; } }
		public override string DescriptionName { get { return "McBryde-Thomas Flat-Polar Sine (No. 1)"; } }

		public override PJ Init()
		{
			return setup(1.48875, 1.36509, false);
		}
	}

	class PJ_fouc : PJ_sts
	{
		public override string Name { get { return "fouc"; } }
		public override string DescriptionName { get { return "Foucaut"; } }

		public override PJ Init()
		{
			return setup(2, 2, true);
		}
	}
}
