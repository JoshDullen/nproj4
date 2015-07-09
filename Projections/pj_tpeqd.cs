using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_tpeqd : PJ
	{
		protected double cp1, sp1, cp2, sp2, ccs, cs, sc, r2z0, z02, dlam2, hz0, thz0, rhshz0, ca, sa, lp, lamc;

		public override string Name { get { return "tpeqd"; } }
		public override string DescriptionName { get { return "Two Point Equidistant"; } }
		public override string DescriptionType { get { return "Misc Sph"; } }
		public override string DescriptionParameters { get { return "lat_1= lon_1= lat_2= lon_2="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +lat_1={0}", Math.Atan2(sp1, cp1)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lon_1={0}", Proj.adjlon(lam0-dlam2)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lat_2={0}", Math.Atan2(sp2, cp2)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +lon_2={0}", Proj.adjlon(lam0+dlam2)*Proj.RAD_TO_DEG);
				return ret.ToString();
			}
		}

		// sphere
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sp=Math.Sin(lp.phi);
			double cp=Math.Cos(lp.phi);
			double dl1=lp.lam+dlam2;
			double dl2=lp.lam-dlam2;
			double z1=Proj.aacos(ctx, sp1*sp+cp1*cp*Math.Cos(dl1));
			double z2=Proj.aacos(ctx, sp2*sp+cp2*cp*Math.Cos(dl2));
			z1*=z1;
			z2*=z2;
			double t=z1-z2;
			xy.x=r2z0*t;
			t=z02-t;
			xy.y=r2z0*Proj.asqrt(4.0*z02*z2-t*t);
			if((ccs*sp-cp*(cs*Math.Sin(dl1)-sc*Math.Sin(dl2)))<0.0) xy.y=-xy.y;

			return xy;
		}

		// sphere
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double cz1=Math.Cos(Libc.hypot(xy.y, xy.x+hz0));
			double cz2=Math.Cos(Libc.hypot(xy.y, xy.x-hz0));
			double s=cz1+cz2;
			double d=cz1-cz2;
			lp.lam=-Math.Atan2(d, (s*thz0));
			lp.phi=Proj.aacos(ctx, Libc.hypot(thz0*s, d)*rhshz0);
			if(xy.y<0.0) lp.phi=-lp.phi;

			// lam--phi now in system relative to P1--P2 base equator
			double sp=Math.Sin(lp.phi);
			double cp=Math.Cos(lp.phi);
			lp.lam-=this.lp;
			s=Math.Cos(lp.lam);
			lp.phi=Proj.aasin(ctx, sa*sp+ca*cp*s);
			lp.lam=Math.Atan2(cp*Math.Sin(lp.lam), sa*cp*s-ca*sp)+lamc;

			return lp;
		}

		public override PJ Init()
		{
			// get control point locations
			double phi_1=Proj.pj_param_r(ctx, parameters, "lat_1");
			double lam_1=Proj.pj_param_r(ctx, parameters, "lon_1");
			double phi_2=Proj.pj_param_r(ctx, parameters, "lat_2");
			double lam_2=Proj.pj_param_r(ctx, parameters, "lon_2");
			if(phi_1==phi_2&&lam_1==lam_2) { Proj.pj_ctx_set_errno(ctx, -25); return null; }
			lam0=Proj.adjlon(0.5*(lam_1+lam_2));
			dlam2=Proj.adjlon(lam_2-lam_1);
			cp1=Math.Cos(phi_1);
			cp2=Math.Cos(phi_2);
			sp1=Math.Sin(phi_1);
			sp2=Math.Sin(phi_2);
			cs=cp1*sp2;
			sc=sp1*cp2;
			ccs=cp1*cp2*Math.Sin(dlam2);
			z02=Proj.aacos(ctx, sp1*sp2+cp1*cp2*Math.Cos(dlam2));
			hz0=0.5*z02;
			double A12=Math.Atan2(cp2*Math.Sin(dlam2), cp1*sp2-sp1*cp2*Math.Cos(dlam2));
			double pp=Proj.aasin(ctx, cp1*Math.Sin(A12));
			ca=Math.Cos(pp);
			sa=Math.Sin(pp);
			lp=Proj.adjlon(Math.Atan2(cp1*Math.Cos(A12), sp1)-hz0);
			dlam2*=0.5;
			lamc=Proj.HALFPI-Math.Atan2(Math.Sin(A12)*sp1, Math.Cos(A12))-dlam2;
			thz0=Math.Tan(hz0);
			rhshz0=0.5/Math.Sin(hz0);
			r2z0=0.5/z02;
			z02*=z02;
			inv=s_inverse;
			fwd=s_forward;
			es=0.0;

			return this;
		}
	}
}
