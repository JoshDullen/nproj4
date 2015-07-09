using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_laea : PJ
	{
		protected double sinb1, cosb1, xmf, ymf, mmf, qp, dd, rq;
		protected double[] apa;
		protected laea_mode mode;

		public override string Name { get { return "laea"; } }
		public override string DescriptionName { get { return "Lambert Azimuthal Equal Area"; } }
		public override string DescriptionType { get { return "Azi, Sph&Ell"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double EPS10=1.0e-10;
		const int NITER=20;
		const double CONV=1.0e-10;

		protected enum laea_mode
		{
			N_POLE=0,
			S_POLE=1,
			EQUIT=2,
			OBLIQ=3
		}

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double sinb=0.0, cosb=0.0, b=0.0;
			double coslam=Math.Cos(lp.lam);
			double sinlam=Math.Sin(lp.lam);
			double sinphi=Math.Sin(lp.phi);
			double q=Proj.pj_qsfn(sinphi, e, one_es);

			if(mode==laea_mode.OBLIQ||mode==laea_mode.EQUIT)
			{
				sinb=q/qp;
				cosb=Math.Sqrt(1.0-sinb*sinb);
			}

			switch(mode)
			{
				case laea_mode.OBLIQ: b=1.0+sinb1*sinb+cosb1*cosb*coslam; break;
				case laea_mode.EQUIT: b=1.0+cosb*coslam; break;
				case laea_mode.N_POLE: b=Proj.HALFPI+lp.phi; q=qp-q; break;
				case laea_mode.S_POLE: b=lp.phi-Proj.HALFPI; q=qp+q; break;
			}

			if(Math.Abs(b)<EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			switch(mode)
			{
				case laea_mode.OBLIQ:
					xy.y=ymf*(b=Math.Sqrt(2.0/b))*(cosb1*sinb-sinb1*cosb*coslam);
					xy.x=xmf*b*cosb*sinlam;
					break;
				case laea_mode.EQUIT:
					xy.y=(b=Math.Sqrt(2.0/(1.0+cosb*coslam)))*sinb*ymf;
					xy.x=xmf*b*cosb*sinlam;
					break;
				case laea_mode.N_POLE:
				case laea_mode.S_POLE:
					if(q>=0.0)
					{
						xy.x=(b=Math.Sqrt(q))*sinlam;
						xy.y=coslam*(mode==laea_mode.S_POLE?b:-b);
					}
					else xy.x=xy.y=0.0;
					break;
			}

			return xy;
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
				case laea_mode.EQUIT:
				case laea_mode.OBLIQ:
					xy.y=1.0+(mode==laea_mode.EQUIT?cosphi*coslam:sinb1*sinphi+cosb1*cosphi*coslam);
					if(xy.y<=EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					xy.x=(xy.y=Math.Sqrt(2.0/xy.y))*cosphi*Math.Sin(lp.lam);
					xy.y*=(mode==laea_mode.EQUIT?sinphi:cosb1*sinphi-sinb1*cosphi*coslam);
					break;
				case laea_mode.N_POLE:
				case laea_mode.S_POLE:
					if(mode==laea_mode.N_POLE) coslam=-coslam;
					if(Math.Abs(lp.phi+phi0)<EPS10) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					xy.y=Proj.FORTPI-lp.phi*0.5;
					xy.y=2.0*(mode==laea_mode.S_POLE?Math.Cos(xy.y):Math.Sin(xy.y));
					xy.x=xy.y*Math.Sin(lp.lam);
					xy.y*=coslam;
					break;
			}

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double ab=0.0;

			switch(mode)
			{
				case laea_mode.EQUIT:
				case laea_mode.OBLIQ:
					xy.x/=dd;
					xy.y*=dd;
					double rho=Libc.hypot(xy.x, xy.y);
					if(rho<EPS10)
					{
						lp.lam=0.0;
						lp.phi=phi0;
						return lp;
					}

					double sCe=2.0*Math.Asin(0.5*rho/rq);
					double cCe=Math.Cos(sCe);
					sCe=Math.Sin(sCe);
					xy.x*=sCe;

					if(mode==laea_mode.OBLIQ)
					{
						ab=cCe*sinb1+xy.y*sCe*cosb1/rho;
						xy.y=rho*cosb1*cCe-xy.y*sinb1*sCe;
					}
					else
					{
						ab=xy.y*sCe/rho;
						xy.y=rho*cCe;
					}
					break;
				case laea_mode.N_POLE:
				case laea_mode.S_POLE:
					if(mode==laea_mode.N_POLE) xy.y=-xy.y;
					double q=xy.x*xy.x+xy.y*xy.y;
					if(q==0)
					{
						lp.lam=0.0;
						lp.phi=phi0;
						return lp;
					}

					//q=P.qp-q;
					ab=1.0-q/qp;
					if(mode==laea_mode.S_POLE) ab=-ab;
					break;
			}

			lp.lam=Math.Atan2(xy.x, xy.y);
			lp.phi=Proj.pj_authlat(Math.Asin(ab), apa);

			return lp;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double cosz=0.0, sinz=0.0;

			double rh=Libc.hypot(xy.x, xy.y);
			lp.phi=rh*0.5;
			if(lp.phi>1.0) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

			lp.phi=2.0*Math.Asin(lp.phi);
			if(mode==laea_mode.OBLIQ||mode==laea_mode.EQUIT)
			{
				sinz=Math.Sin(lp.phi);
				cosz=Math.Cos(lp.phi);
			}

			switch(mode)
			{
				case laea_mode.EQUIT:
					lp.phi=Math.Abs(rh)<=EPS10?0.0:Math.Asin(xy.y*sinz/rh);
					xy.x*=sinz;
					xy.y=cosz*rh;
					break;
				case laea_mode.OBLIQ:
					lp.phi=Math.Abs(rh)<=EPS10?phi0:Math.Asin(cosz*sinb1+xy.y*sinz*cosb1/rh);
					xy.x*=sinz*cosb1;
					xy.y=(cosz-Math.Sin(lp.phi)*sinb1)*rh;
					break;
				case laea_mode.N_POLE:
					xy.y=-xy.y;
					lp.phi=Proj.HALFPI-lp.phi;
					break;
				case laea_mode.S_POLE:
					lp.phi-=Proj.HALFPI;
					break;
			}

			lp.lam=(xy.y==0.0&&(mode==laea_mode.EQUIT||mode==laea_mode.OBLIQ))?0.0:Math.Atan2(xy.x, xy.y);

			return lp;
		}

		public override PJ Init()
		{
			double t=Math.Abs(phi0);

			if(Math.Abs(t-Proj.HALFPI)<EPS10) mode=phi0<0.0?laea_mode.S_POLE:laea_mode.N_POLE;
			else if(Math.Abs(t)<EPS10) mode=laea_mode.EQUIT;
			else mode=laea_mode.OBLIQ;

			if(es!=0)
			{
				e=Math.Sqrt(es);
				qp=Proj.pj_qsfn(1.0, e, one_es);
				mmf=0.5/(1.0-es);
				apa=Proj.pj_authset(es);
				switch(mode)
				{
					case laea_mode.N_POLE:
					case laea_mode.S_POLE:
						dd=1.0;
						break;
					case laea_mode.EQUIT:
						rq=Math.Sqrt(0.5*qp);
						dd=1.0/rq;
						xmf=1.0;
						ymf=0.5*qp;
						break;
					case laea_mode.OBLIQ:
						rq=Math.Sqrt(0.5*qp);
						double sinphi=Math.Sin(phi0);
						sinb1=Proj.pj_qsfn(sinphi, e, one_es)/qp;
						cosb1=Math.Sqrt(1.0-sinb1*sinb1);
						dd=Math.Cos(phi0)/(Math.Sqrt(1.0-es*sinphi*sinphi)*rq*cosb1);
						xmf=rq;
						ymf=xmf/dd;
						xmf*=dd;
						break;
				}
				inv=e_inverse;
				fwd=e_forward;
			}
			else
			{
				if(mode==laea_mode.OBLIQ)
				{
					sinb1=Math.Sin(phi0);
					cosb1=Math.Cos(phi0);
				}
				inv=s_inverse;
				fwd=s_forward;
			}

			return this;
		}
	}
}
