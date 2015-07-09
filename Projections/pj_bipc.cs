using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Projections
{
	class PJ_bipc : PJ
	{
		protected bool noskew;

		public override string Name { get { return "bipc"; } }
		public override string DescriptionName { get { return "Bipolar conic of western hemisphere"; } }
		public override string DescriptionType { get { return "Conic, Sph"; } }
		public override string DescriptionParameters { get { return "ns"; } }
		public override bool Invertible { get { return true; } }

		protected override string Proj4ParameterString
		{
			get
			{
				return noskew?" +ns":"";
			}
		}

		const double EPS=1.0e-10;
		const double EPS10=1.0e-10;
		const double ONEEPS=1.000000001;
		const int NITER=10;
		const double lamB=-0.34894976726250681539;
		const double n=0.63055844881274687180;
		const double F=1.89724742567461030582;
		const double Azab=0.81650043674686363166;
		const double Azba=1.82261843856185925133;
		const double T=1.27246578267089012270;
		const double rhoc=1.20709121521568721927;
		const double cAzc=0.69691523038678375519;
		const double sAzc=0.71715351331143607555;
		const double C45=0.70710678118654752469;
		const double S45=0.70710678118654752410;
		const double C20=0.93969262078590838411;
		const double S20=-0.34202014332566873287;
		const double R110=1.91986217719376253360;
		const double R104=1.81514242207410275904;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double cphi=Math.Cos(lp.phi);
			double sphi=Math.Sin(lp.phi);
			double sdlam=lamB-lp.lam;
			double cdlam=Math.Cos(sdlam);
			sdlam=Math.Sin(sdlam);

			double tphi, Az;
			if(Math.Abs(Math.Abs(lp.phi)-Proj.HALFPI)<EPS10)
			{
				Az=lp.phi<0.0?Proj.PI:0.0;
				tphi=Libc.HUGE_VAL;
			}
			else
			{
				tphi=sphi/cphi;
				Az=Math.Atan2(sdlam, C45*(tphi-cdlam));
			}

			bool tag=Az>Azba;

			double z, Av;

			if(tag)
			{
				cdlam=Math.Cos(sdlam=lp.lam+R110);
				sdlam=Math.Sin(sdlam);
				z=S20*sphi+C20*cphi*cdlam;
				if(Math.Abs(z)>1.0)
				{
					if(Math.Abs(z)>ONEEPS) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					else z=z<0.0?-1.0:1.0;
				}
				else z=Math.Acos(z);

				if(tphi!=Libc.HUGE_VAL) Az=Math.Atan2(sdlam, (C20*tphi-S20*cdlam));
				Av=Azab;
				xy.y=rhoc;
			}
			else
			{
				z=S45*(sphi+cphi*cdlam);
				if(Math.Abs(z)>1.0)
				{
					if(Math.Abs(z)>ONEEPS) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
					else z=z<0.0?-1.0:1.0;
				}
				else z=Math.Acos(z);

				Av=Azba;
				xy.y=-rhoc;
			}

			if(z<0.0) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			double t=Math.Pow(Math.Tan(0.5*z), n);
			double r=F*t;
			double al=0.5*(R104-z);

			if(al<0.0) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			al=(t+Math.Pow(al, n))/T;
			if(Math.Abs(al)>1.0)
			{
				if(Math.Abs(al)>ONEEPS) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
				else al=al<0.0?-1.0:1.0;
			}
			else al=Math.Acos(al);

			t=n*(Av-Az);
			if(Math.Abs(t)<al) r/=Math.Cos(al+(tag?t:-t));
			xy.x=r*Math.Sin(t);
			xy.y+=(tag?-r:r)*Math.Cos(t);

			if(noskew)
			{
				t=xy.x;
				xy.x=-xy.x*cAzc-xy.y*sAzc;
				xy.y=-xy.y*cAzc+t*sAzc;
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			if(noskew)
			{
				double t=xy.x;
				xy.x=-xy.x*cAzc+xy.y*sAzc;
				xy.y=-xy.y*cAzc-t*sAzc;
			}

			double s, c, Av;

			bool neg=(xy.x<0.0);
			if(neg)
			{
				xy.y=rhoc-xy.y;
				s=S20;
				c=C20;
				Av=Azab;
			}
			else
			{
				xy.y+=rhoc;
				s=S45;
				c=C45;
				Av=Azba;
			}

			double rl=Libc.hypot(xy.x, xy.y);
			double rp=rl, r=rl;
			double Az=Math.Atan2(xy.x, xy.y);
			double fAz=Math.Abs(Az);

			double z=0;
			int i=NITER;
			for(; i>0; i--)
			{
				z=2.0*Math.Atan(Math.Pow(r/F, 1/n));
				double al=Math.Acos((Math.Pow(Math.Tan(0.5*z), n)+Math.Pow(Math.Tan(0.5*(R104-z)), n))/T);
				if(fAz<al) r=rp*Math.Cos(al+(neg?Az:-Az));
				if(Math.Abs(rl-r)<EPS) break;
				rl=r;
			}

			if(i==0) { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

			Az=Av-Az/n;
			lp.phi=Math.Asin(s*Math.Cos(z)+c*Math.Sin(z)*Math.Cos(Az));
			lp.lam=Math.Atan2(Math.Sin(Az), c/Math.Tan(z)-s*Math.Cos(Az));

			if(neg) lp.lam-=R110;
			else lp.lam=lamB-lp.lam;

			return lp;
		}

		public override PJ Init()
		{
			noskew=Proj.pj_param_b(ctx, parameters, "ns");
			inv=s_inverse;
			fwd=s_forward;
			es=0.0;

			return this;
		}
	}
}
