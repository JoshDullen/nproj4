using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_urmfps : PJ
	{
		protected double n, C_y;

		public override string Name { get { return "urmfps"; } }
		public override string DescriptionName { get { return "Urmaev Flat-Polar Sinusoidal"; } }
		public override string DescriptionType { get { return "PCyl, Sph"; } }
		public override string DescriptionParameters { get { return "n="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +n={0}", n);
				return ret.ToString();
			}
		}

		const double C_x=0.8773826753;
		const double Cy=1.139753528477;

		// sphere
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			lp.phi=Proj.aasin(ctx, n*Math.Sin(lp.phi));
			xy.x=C_x*lp.lam*Math.Cos(lp.phi);
			xy.y=C_y*lp.phi;
			return xy;
		}

		// sphere
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y/=C_y;
			lp.phi=Proj.aasin(ctx, Math.Sin(xy.y)/n);
			lp.lam=xy.x/(C_x*Math.Cos(xy.y));

			return lp;
		}

		protected PJ_urmfps setup()
		{
			C_y=Cy/n;
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}

		public override PJ Init()
		{
			if(Proj.pj_param_t(ctx, parameters, "n"))
			{
				n=Proj.pj_param_d(ctx, parameters, "n");
				if(n<=0.0||n>1.0) { Proj.pj_ctx_set_errno(ctx, -40); return null; }
			}
			else { Proj.pj_ctx_set_errno(ctx, -40); return null; }

			return setup();
		}
	}

	class PJ_wag1 : PJ_urmfps
	{
		public override string Name { get { return "wag1"; } }
		public override string DescriptionName { get { return "Wagner I (Kavraisky VI)"; } }
		public override string DescriptionParameters { get { return ""; } }

		protected override string Proj4ParameterString
		{
			get
			{
				return "";
			}
		}

		public override PJ Init()
		{
			n=0.8660254037844386467637231707;

			return setup();
		}
	}
}
