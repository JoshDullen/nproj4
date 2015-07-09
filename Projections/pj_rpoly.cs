using System;
using System.Text;

namespace Free.Ports.Proj4
{
	class PJ_rpoly : PJ
	{
		protected double phi1, fxa, fxb;
		protected bool mode;

		public override string Name { get { return "rpoly"; } }
		public override string DescriptionName { get { return "Rectangular Polyconic"; } }
		public override string DescriptionType { get { return "Conic, Sph, no inv."; } }
		public override string DescriptionParameters { get { return "lat_ts="; } }
		public override bool Invertible { get { return false; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_ts={0}", phi1*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		const double EPS=1.0e-9;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;
			double fa;

			if(mode) fa=Math.Tan(lp.lam*fxb)*fxa;
			else fa=0.5*lp.lam;

			if(Math.Abs(lp.phi)<EPS)
			{
				xy.x=fa+fa;
				xy.y=-phi0;
			}
			else
			{
				xy.y=1.0/Math.Tan(lp.phi);
				fa=2.0*Math.Atan(fa*Math.Sin(lp.phi));
				xy.x=Math.Sin(fa)*xy.y;
				xy.y=lp.phi-phi0+(1.0-Math.Cos(fa))*xy.y;
			}

			return xy;
		}

		public override PJ Init()
		{
			phi1=Math.Abs(Proj.pj_param_r(ctx, parameters, "lat_ts"));
			mode=phi1>EPS;

			if(mode)
			{
				fxb=0.5*Math.Sin(phi1);
				fxa=0.5/fxb;
			}

			es=0.0;
			fwd=s_forward;

			return this;
		}
	}
}
