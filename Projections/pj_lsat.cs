using System;
using System.Text;
using Free.Ports.Proj4.LibcStuff;

// based upon Snyder and Linck, USGS-NMD
namespace Free.Ports.Proj4.Projections
{
	class PJ_lsat : PJ
	{
		protected double a2, a4, b, c1, c3, q, t, u, w, p22, sa, ca, xj, rlm, rlm2;

		public override string Name { get { return "lsat"; } }
		public override string DescriptionName { get { return "Space oblique for LANDSAT"; } }
		public override string DescriptionType { get { return "Cyl, Sph&Ell"; } }
		public override string DescriptionParameters { get { return "lsat= path="; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				StringBuilder ret=new StringBuilder();
				double alf=Math.Asin(sa);
				if(Math.Abs(p22*1440.0-103.2669323)<1.0e-10&&Math.Abs(alf-Proj.DEG_TO_RAD*99.092)<1.0e-10)
				{
					ret.Append(" +lsat=1");
					int path=(int)Math.Floor((lam0+Proj.DEG_TO_RAD*128.87)/(7-Proj.TWOPI/251.0)+0.5);
					ret.AppendFormat(nc, " +path={0}", path);
				}
				else if(Math.Abs(p22*1440.0-98.8841202)<1.0e-10&&Math.Abs(alf-Proj.DEG_TO_RAD*98.2)<1.0e-10)
				{
					ret.Append(" +lsat=4");
					int path=(int)Math.Floor((lam0+Proj.DEG_TO_RAD*129.3)/(7-Proj.TWOPI/233.0)+0.5);
					ret.AppendFormat(nc, " +path={0}", path);
				}
				else throw new InvalidOperationException("Project settings invalid, or incomplete.");

				return ret.ToString();
			}
		}

		const double TOL=1.0e-7;
		const double PI_HALFPI=4.71238898038468985766;
		const double TWOPI_HALFPI=7.85398163397448309610;

		void seraz0(double lam, double mult)
		{
			lam*=Proj.DEG_TO_RAD;
			double sd=Math.Sin(lam);
			double sdsq=sd*sd;
			double s=p22*sa*Math.Cos(lam)*Math.Sqrt((1.0+t*sdsq)/((1.0+w*sdsq)*(1.0+q*sdsq)));
			double d__1=1.0+q*sdsq;
			double h=Math.Sqrt((1.0+q*sdsq)/(1.0+w*sdsq))*((1.0+w*sdsq)/(d__1*d__1)-p22*ca);
			double sq=Math.Sqrt(xj*xj+s*s);
			double fc=mult*(h*xj-s*s)/sq;
			b+=fc;
			a2+=fc*Math.Cos(lam+lam);
			a4+=fc*Math.Cos(lam*4.0);
			fc=mult*s*(h+xj)/sq;
			c1+=fc*Math.Cos(lam);
			c3+=fc*Math.Cos(lam*3.0);
		}

		// ellipsoid
		XY e_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			if(lp.phi>Proj.HALFPI) lp.phi=Proj.HALFPI;
			else if(lp.phi<-Proj.HALFPI) lp.phi=-Proj.HALFPI;
			double lampp=lp.phi>=0.0?Proj.HALFPI:PI_HALFPI;
			double tanphi=Math.Tan(lp.phi);

			double lamt=0, lamdp=0;
			int l;
			for(int nn=0; ; )
			{
				double sav=lampp;
				double lamtp=lp.lam+p22*lampp;
				double cl=Math.Cos(lamtp);
				if(Math.Abs(cl)<TOL) lamtp-=TOL;
				double fac=lampp-Math.Sin(lampp)*(cl<0.0?-Proj.HALFPI:Proj.HALFPI);
				for(l=50; l>0; l--)
				{
					lamt=lp.lam+p22*sav;
					double c=Math.Cos(lamt);
					if(Math.Abs(c)<TOL) lamt-=TOL;
					double xlam=(one_es*tanphi*sa+Math.Sin(lamt)*ca)/c;
					lamdp=Math.Atan(xlam)+fac;
					if(Math.Abs(Math.Abs(sav)-Math.Abs(lamdp))<TOL) break;
					sav=lamdp;
				}
				nn++;
				if(l==0||nn>=3||(lamdp>rlm&&lamdp<rlm2)) break;
				if(lamdp<=rlm) lampp=TWOPI_HALFPI;
				else if(lamdp>=rlm2) lampp=Proj.HALFPI;
			}

			if(l!=0)
			{
				double sp=Math.Sin(lp.phi);
				double phidp=Proj.aasin(ctx, (one_es*ca*sp-sa*Math.Cos(lp.phi)*Math.Sin(lamt))/Math.Sqrt(1.0-es*sp*sp));
				double tanph=Math.Log(Math.Tan(Proj.FORTPI+0.5*phidp));
				double sd=Math.Sin(lamdp);
				double sdsq=sd*sd;
				double s=p22*sa*Math.Cos(lamdp)*Math.Sqrt((1.0+t*sdsq)/((1.0+w*sdsq)*(1.0+q*sdsq)));
				double d=Math.Sqrt(xj*xj+s*s);
				xy.x=b*lamdp+a2*Math.Sin(2.0*lamdp)+a4*Math.Sin(lamdp*4.0)-tanph*s/d;
				xy.y=c1*sd+c3*Math.Sin(lamdp*3.0)+tanph*xj/d;
			}
			else xy.x=xy.y=Libc.HUGE_VAL;

			return xy;
		}

		// ellipsoid
		LP e_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double s;

			double lamdp=xy.x/b;
			int nn=50;
			double sav;
			do
			{
				sav=lamdp;
				double sd=Math.Sin(lamdp);
				double sdsq=sd*sd;
				s=p22*sa*Math.Cos(lamdp)*Math.Sqrt((1.0+t*sdsq)/((1.0+w*sdsq)*(1.0+q*sdsq)));
				lamdp=xy.x+xy.y*s/xj-a2*Math.Sin(2.0*lamdp)-a4*Math.Sin(lamdp*4.0)-s/xj*(c1*Math.Sin(lamdp)+c3*Math.Sin(lamdp*3.0));
				lamdp/=b;
				nn--;
			} while(Math.Abs(lamdp-sav)>=TOL&&nn>0);

			double sl=Math.Sin(lamdp);
			double fac=Math.Exp(Math.Sqrt(1.0+s*s/xj/xj)*(xy.y-c1*sl-c3*Math.Sin(lamdp*3.0)));
			double phidp=2.0*(Math.Atan(fac)-Proj.FORTPI);
			double dd=sl*sl;
			if(Math.Abs(Math.Cos(lamdp))<TOL) lamdp-=TOL;
			double spp=Math.Sin(phidp);
			double sppsq=spp*spp;
			double lamt=Math.Atan(((1.0-sppsq*rone_es)*Math.Tan(lamdp)*ca-spp*sa*Math.Sqrt((1.0+q*dd)*(1.0-sppsq)-sppsq*u)/Math.Cos(lamdp))/(1.0-sppsq*(1.0+u)));
			sl=lamt>=0.0?1.0:-1.0;
			double scl=Math.Cos(lamdp)>=0.0?1.0:-1;
			lamt-=Proj.HALFPI*(1.0-scl)*sl;
			lp.lam=lamt-p22*lamdp;
			if(Math.Abs(sa)<TOL) lp.phi=Proj.aasin(ctx, spp/Math.Sqrt(one_es*one_es+es*sppsq));
			else lp.phi=Math.Atan((Math.Tan(lamdp)*Math.Cos(lamt)-ca*Math.Sin(lamt))/(one_es*sa));

			return lp;
		}

		public override PJ Init()
		{
			int land=Proj.pj_param_i(ctx, parameters, "lsat");
			if(land<=0||land>5) { Proj.pj_ctx_set_errno(ctx, -28); return null; }

			int path=Proj.pj_param_i(ctx, parameters, "path");
			if(path<=0||path>(land<=3?251:233)) { Proj.pj_ctx_set_errno(ctx, -29); return null; }
			double alf;
			if(land<=3)
			{
				lam0=Proj.DEG_TO_RAD*128.87-Proj.TWOPI/251.0*path;
				p22=103.2669323;
				alf=Proj.DEG_TO_RAD*99.092;
			}
			else
			{
				lam0=Proj.DEG_TO_RAD*129.3-Proj.TWOPI/233.0*path;
				p22=98.8841202;
				alf=Proj.DEG_TO_RAD*98.2;
			}

			p22/=1440.0;
			sa=Math.Sin(alf);
			ca=Math.Cos(alf);
			if(Math.Abs(ca)<1e-9) ca=1e-9;
			double esc=es*ca*ca;
			double ess=es*sa*sa;
			w=(1.0-esc)*rone_es;
			w=w*w-1.0;
			q=ess*rone_es;
			t=ess*(2.0-es)*rone_es*rone_es;
			u=esc*rone_es;
			xj=one_es*one_es*one_es;
			rlm=Proj.PI*(1.0/248.0+0.5161290322580645);
			rlm2=rlm+Proj.TWOPI;
			a2=a4=b=c1=c3=0.0;
			seraz0(0.0, 1.0);
			for(double lam=9.0; lam<=81.0001; lam+=18.0) seraz0(lam, 4.0);
			for(double lam=18; lam<=72.0001; lam+=18.0) seraz0(lam, 2.0);
			seraz0(90.0, 1.0);
			a2/=30.0;
			a4/=60.0;
			b/=30.0;
			c1/=15.0;
			c3/=45.0;

			inv=e_inverse;
			fwd=e_forward;

			return this;
		}
	}
}
