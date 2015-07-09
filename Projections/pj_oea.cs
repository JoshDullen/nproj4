using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_oea : PJ
	{
		protected double theta, m, n, two_r_m, two_r_n, rm, rn, hm, hn, cp0, sp0;

		public override string Name { get { return "oea"; } }
		public override string DescriptionName { get { return "Oblated Equal Area"; } }
		public override string DescriptionType { get { return "Misc Sph"; } }
		public override string DescriptionParameters { get { return "n= m= theta="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +n={0}", n);
				ret.AppendFormat(nc, " +m={0}", m);
				ret.AppendFormat(nc, " +theta={0}", theta*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		// sphere
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double cp=Math.Cos(lp.phi);
			double sp=Math.Sin(lp.phi);
			double cl=Math.Cos(lp.lam);
			double Az=Proj.aatan2(cp*Math.Sin(lp.lam), cp0*sp-sp0*cp*cl)+theta;
			double shz=Math.Sin(0.5*Proj.aacos(ctx, sp0*sp+cp0*cp*cl));
			double M=Proj.aasin(ctx, shz*Math.Sin(Az));
			double N=Proj.aasin(ctx, shz*Math.Cos(Az)*Math.Cos(M)/Math.Cos(M*two_r_m));
			xy.y=n*Math.Sin(N*two_r_n);
			xy.x=m*Math.Sin(M*two_r_m)*Math.Cos(N)/Math.Cos(N*two_r_n);

			return xy;
		}

		// sphere
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double N=hn*Proj.aasin(ctx, xy.y*rn);
			double M=hm*Proj.aasin(ctx, xy.x*rm*Math.Cos(N*two_r_n)/Math.Cos(N));
			double xp=2.0*Math.Sin(M);
			double yp=2.0*Math.Sin(N)*Math.Cos(M*two_r_m)/Math.Cos(M);
			double Az=Proj.aatan2(xp, yp)-theta;
			double cAz=Math.Cos(Az);
			double z=2.0*Proj.aasin(ctx, 0.5*Libc.hypot(xp, yp));
			double sz=Math.Sin(z);
			double cz=Math.Cos(z);
			lp.phi=Proj.aasin(ctx, sp0*cz+cp0*sz*cAz);
			lp.lam=Proj.aatan2(sz*Math.Sin(Az), cp0*cz-sp0*sz*cAz);

			return lp;
		}

		public override PJ Init()
		{
			n=Proj.pj_param_d(ctx, parameters, "n");
			m=Proj.pj_param_d(ctx, parameters, "m");

			if(n<=0.0||m<=0.0) { Proj.pj_ctx_set_errno(ctx, -39); return null; }

			theta=Proj.pj_param_r(ctx, parameters, "theta");
			sp0=Math.Sin(phi0);
			cp0=Math.Cos(phi0);
			rn=1.0/n;
			rm=1.0/m;
			two_r_n=2.0*rn;
			two_r_m=2.0*rm;
			hm=0.5*m;
			hn=0.5*n;
			fwd=s_forward;
			inv=s_inverse;
			es=0.0;

			return this;
		}
	}
}
