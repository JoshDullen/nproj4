using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_nsper : PJ
	{
		protected double height, sinph0, cosph0, p, rp, pn1, pfact, h, cg, sg, sw, cw;
		protected nsper_mode mode;
		protected bool tilt;

		public override string Name { get { return "nsper"; } }
		public override string DescriptionName { get { return "Near-sided perspective"; } }
		public override string DescriptionType { get { return "Azi, Sph"; } }
		public override string DescriptionParameters { get { return "h="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +h={0}", height);
				return ret.ToString();
			}
		}

		const double EPS10=1.0e-10;

		protected enum nsper_mode
		{
			N_POLE=0,
			S_POLE=1,
			EQUIT=2,
			OBLIQ=3
		}

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sinphi=Math.Sin(lp.phi);
			double cosphi=Math.Cos(lp.phi);
			double coslam=Math.Cos(lp.lam);

			switch(mode)
			{
				case nsper_mode.OBLIQ: xy.y=sinph0*sinphi+cosph0*cosphi*coslam; break;
				case nsper_mode.EQUIT: xy.y=cosphi*coslam; break;
				case nsper_mode.S_POLE: xy.y=-sinphi; break;
				case nsper_mode.N_POLE: xy.y=sinphi; break;
			}

			if(xy.y<rp) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			xy.y=pn1/(p-xy.y);
			xy.x=xy.y*cosphi*Math.Sin(lp.lam);

			switch(mode)
			{
				case nsper_mode.OBLIQ: xy.y*=(cosph0*sinphi-sinph0*cosphi*coslam); break;
				case nsper_mode.EQUIT: xy.y*=sinphi; break;
				case nsper_mode.N_POLE: coslam=-coslam; xy.y*=cosphi*coslam; break;
				case nsper_mode.S_POLE: xy.y*=cosphi*coslam; break;
			}

			if(tilt)
			{
				double yt=xy.y*cg+xy.x*sg;
				double ba=1.0/(yt*sw*h+cw);
				xy.x=(xy.x*cg-xy.y*sg)*cw*ba;
				xy.y=yt*ba;
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			if(tilt)
			{
				double yt=1.0/(pn1-xy.y*sw);
				double bm=pn1*xy.x*yt;
				double bq=pn1*xy.y*cw*yt;
				xy.x=bm*cg+bq*sg;
				xy.y=bq*cg-bm*sg;
			}

			double rh=Libc.hypot(xy.x, xy.y);
			double sinz=1.0-rh*rh*pfact;
			if(sinz<0.0) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

			sinz=(p-Math.Sqrt(sinz))/(pn1/rh+rh/pn1);
			double cosz=Math.Sqrt(1.0-sinz*sinz);

			if(Math.Abs(rh)<=EPS10)
			{
				lp.lam=0.0;
				lp.phi=phi0;
			}
			else
			{
				switch(mode)
				{
					case nsper_mode.OBLIQ:
						lp.phi=Math.Asin(cosz*sinph0+xy.y*sinz*cosph0/rh);
						xy.y=(cosz-sinph0*Math.Sin(lp.phi))*rh;
						xy.x*=sinz*cosph0;
						break;
					case nsper_mode.EQUIT: lp.phi=Math.Asin(xy.y*sinz/rh); xy.y=cosz*rh; xy.x*=sinz; break;
					case nsper_mode.N_POLE: lp.phi=Math.Asin(cosz); xy.y=-xy.y; break;
					case nsper_mode.S_POLE: lp.phi=-Math.Asin(cosz); break;
				}

				lp.lam=Math.Atan2(xy.x, xy.y);
			}

			return lp;
		}

		protected PJ_nsper setup()
		{
			height=Proj.pj_param_d(ctx, parameters, "h");
			if(height<=0.0) { Proj.pj_ctx_set_errno(ctx, -30); return null; }
			if(Math.Abs(Math.Abs(phi0)-Proj.HALFPI)<EPS10) mode=phi0<0.0?nsper_mode.S_POLE:nsper_mode.N_POLE;
			else if(Math.Abs(phi0)<EPS10) mode=nsper_mode.EQUIT;
			else
			{
				mode=nsper_mode.OBLIQ;
				sinph0=Math.Sin(phi0);
				cosph0=Math.Cos(phi0);
			}
			pn1=height/a; // normalize by radius
			p=1.0+pn1;
			rp=1.0/p;
			h=1.0/pn1;
			pfact=(p+1.0)*h;
			inv=s_inverse;
			fwd=s_forward;
			es=0.0;

			return this;
		}

		public override PJ Init()
		{
			tilt=false;

			return setup();
		}
	}

	class PJ_tpers : PJ_nsper
	{
		public override string Name { get { return "tpers"; } }
		public override string DescriptionName { get { return "Tilted perspective"; } }
		public override string DescriptionParameters { get { return "tilt= azi= h="; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				ret.AppendFormat(nc, " +tilt={0}", Math.Atan2(sg, cg)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +azi={0}", Math.Atan2(sw, cw)*Proj.RAD_TO_DEG);
				ret.AppendFormat(nc, " +h={0}", height);
				return ret.ToString();
			}
		}

		public override PJ Init()
		{
			double omega=Proj.pj_param_d(ctx, parameters, "tilt")*Proj.DEG_TO_RAD;
			double gamma=Proj.pj_param_d(ctx, parameters, "azi")*Proj.DEG_TO_RAD;
			tilt=true;
			cg=Math.Cos(gamma); sg=Math.Sin(gamma);
			cw=Math.Cos(omega); sw=Math.Sin(omega);

			return setup();
		}
	}
}
