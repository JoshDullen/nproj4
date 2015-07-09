using System;
using System.Text;

namespace Free.Ports.Proj4.Projections
{
	class PJ_lagrng : PJ
	{
		protected double hrw, rw, a1;

		public override string Name { get { return "lagrng"; } }
		public override string DescriptionName { get { return "Lagrange"; } }
		public override string DescriptionType { get { return "Misc Sph, no inv."; } }
		public override string DescriptionParameters { get { return "W= lat_1="; } }
		public override bool Invertible { get { return false; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +W={0}", 1/rw);
				double tmp=Math.Pow(a1, 1/hrw);
				tmp=(1-tmp)/(tmp+1);
				if(tmp>1) tmp=1;
				else if(tmp<-1) tmp=-1;
				ret.AppendFormat(nc, " +lat_1={0}", Math.Asin(tmp)*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		public const double TOL=1.0e-10;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<TOL)
			{
				xy.x=0;
				xy.y=lp.phi<0?-2.0:2.0;
			}
			else
			{
				lp.phi=Math.Sin(lp.phi);
				double v=a1*Math.Pow((1.0+lp.phi)/(1.0-lp.phi), hrw);
				lp.lam*=rw;
				double c=0.5*(v+1.0/v)+Math.Cos(lp.lam);
				if(c<TOL) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
				xy.x=2.0*Math.Sin(lp.lam)/c;
				xy.y=(v-1.0/v)/c;
			}

			return xy;
		}

		public override PJ Init()
		{
			double phi1;

			rw=Proj.pj_param_d(ctx, parameters, "W");
			if(rw<=0) { Proj.pj_ctx_set_errno(ctx, -27); return null; }

			rw=1.0/rw;
			hrw=0.5*rw;
			phi1=Proj.pj_param_r(ctx, parameters, "lat_1");
			phi1=Math.Sin(phi1);
			if(Math.Abs(Math.Abs(phi1)-1.0)<TOL) { Proj.pj_ctx_set_errno(ctx, -22); return null; }
			a1=Math.Pow((1.0-phi1)/(1.0+phi1), hrw);
			es=0.0;
			fwd=s_forward;

			return this;
		}
	}
}
