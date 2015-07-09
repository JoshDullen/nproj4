using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_bonne : PJ
	{
		protected double phi1;
		protected double cphi1;
		protected double am1;
		protected double m1;
		protected double[] en;

		public override string Name { get { return "bonne"; } }
		public override string DescriptionName { get { return "Bonne (Werner lat_1=90)"; } }
		public override string DescriptionType { get { return "Conic, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "lat_1="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", phi1*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		const double EPS10=1.0e-10;

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double E=Math.Sin(lp.phi), c=Math.Cos(lp.phi);

			double rh=am1+m1-Proj.pj_mlfn(lp.phi, E, c, en);
			E=c*lp.lam/(rh*Math.Sqrt(1.0-es*E*E));
			xy.x=rh*Math.Sin(E);
			xy.y=am1-rh*Math.Cos(E);

			return xy;
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double rh=cphi1+phi1-lp.phi;
			if(Math.Abs(rh)>EPS10)
			{
				double E=lp.lam*Math.Cos(lp.phi)/rh;
				xy.x=rh*Math.Sin(E);
				xy.y=cphi1-rh*Math.Cos(E);
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y=cphi1-xy.y;
			double rh=Libc.hypot(xy.x, xy.y);

			lp.phi=cphi1+phi1-rh;
			if(Math.Abs(lp.phi)>Proj.HALFPI) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }
			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<=EPS10) lp.lam=0.0;
			else lp.lam=rh*Math.Atan2(xy.x, xy.y)/Math.Cos(lp.phi);

			return lp;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			xy.y=am1-xy.y;
			double rh=Libc.hypot(xy.x, xy.y);

			lp.phi=Proj.pj_inv_mlfn(ctx, am1+m1-rh, es, en);
			double s=Math.Abs(lp.phi);
			if(s<Proj.HALFPI)
			{
				s=Math.Sin(lp.phi);
				lp.lam=rh*Math.Atan2(xy.x, xy.y)*Math.Sqrt(1.0-es*s*s)/Math.Cos(lp.phi);
			}
			else if(Math.Abs(s-Proj.HALFPI)<=EPS10) lp.lam=0.0;
			else Proj.pj_ctx_set_errno(ctx, -20);

			return lp;
		}

		public override PJ Init()
		{
			phi1=Proj.pj_param_r(ctx, parameters, "lat_1");
			if(Math.Abs(phi1)<EPS10) { Proj.pj_ctx_set_errno(ctx, -23); return null; }
			if(es!=0)
			{
				en=Proj.pj_enfn(es);
				am1=Math.Sin(phi1);
				double c=Math.Cos(phi1);
				m1=Proj.pj_mlfn(phi1, am1, c, en);
				am1=c/(Math.Sqrt(1.0-es*am1*am1)*am1);
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{
				if(Math.Abs(phi1)+EPS10>=Proj.HALFPI) cphi1=0.0;
				else cphi1=1.0/Math.Tan(phi1);
				inv=s_inverse;
				fwd=s_forward;
			}

			return this;
		}
	}
}
