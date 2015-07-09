using System;

namespace Free.Ports.Proj4.Projections
{
	class PJ_vandg : PJ
	{
		public override string Name { get { return "vandg"; } }
		public override string DescriptionName { get { return "van der Grinten (I)"; } }
		public override string DescriptionType { get { return "Misc Sph"; } }
		public override string DescriptionParameters { get { return ""; } }
		public override bool Invertible { get { return true; } }

		const double TOL=1.0e-10;
		const double THIRD=0.33333333333333333333;
		const double C2_27=0.07407407407407407407;
		const double PI4_3=4.18879020478639098458;
		const double PISQ=9.86960440108935861869;
		const double TPISQ=19.73920880217871723738;
		const double HPISQ=4.93480220054467930934;

		// spheroid
		XY s_forward(LP lp)
		{
			XY xy;
			xy.x=xy.y=0;

			double p2=Math.Abs(lp.phi/Proj.HALFPI);
			if((p2-TOL)>1.0) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }

			if(p2>1.0) p2=1.0;

			if(Math.Abs(lp.phi)<=TOL)
			{
				xy.x=lp.lam;
				xy.y=0.0;
			}
			else if(Math.Abs(lp.lam)<=TOL||Math.Abs(p2-1.0)<TOL)
			{
				xy.x=0.0;
				xy.y=Proj.PI*Math.Tan(0.5*Math.Asin(p2));
				if(lp.phi<0.0) xy.y=-xy.y;
			}
			else
			{
				double al=0.5*Math.Abs(Proj.PI/lp.lam-lp.lam/Proj.PI);
				double al2=al*al;
				double g=Math.Sqrt(1.0-p2*p2);
				g=g/(p2+g-1.0);
				double g2=g*g;
				p2=g*(2.0/p2-1.0);
				p2=p2*p2;
				xy.x=g-p2;
				g=p2+al2;
				xy.x=Proj.PI*(al*xy.x+Math.Sqrt(al2*xy.x*xy.x-g*(g2-p2)))/g;
				if(lp.lam<0.0) xy.x=-xy.x;
				xy.y=Math.Abs(xy.x/Proj.PI);
				xy.y=1.0-xy.y*(xy.y+2.0*al);
				if(xy.y<-TOL) { Proj.pj_ctx_set_errno(ctx, -20); return xy; }
				if(xy.y<0.0) xy.y=0.0;
				else xy.y=Math.Sqrt(xy.y)*(lp.phi<0.0?-Proj.PI:Proj.PI);
			}

			return xy;
		}

		// spheroid
		LP s_inverse(XY xy)
		{
			LP lp;
			lp.lam=lp.phi=0;

			double x2=xy.x*xy.x;
			double ay=Math.Abs(xy.y);
			if(ay<TOL)
			{
				lp.phi=0.0;
				double t1=x2*x2+TPISQ*(x2+HPISQ);
				lp.lam=Math.Abs(xy.x)<=TOL?0.0:0.5*(x2-PISQ+Math.Sqrt(t1))/xy.x;
				return lp;
			}

			double y2=xy.y*xy.y;
			double r=x2+y2;
			double r2=r*r;
			double c1=-Proj.PI*ay*(r+PISQ);
			double c3=r2+Proj.TWOPI*(ay*r+Proj.PI*(y2+Proj.PI*(ay+Proj.HALFPI)));
			double c2=c1+PISQ*(r-3.0*y2);
			double c0=Proj.PI*ay;
			c2/=c3;
			double al=c1/c3-THIRD*c2*c2;
			double m=2.0*Math.Sqrt(-THIRD*al);
			double d=C2_27*c2*c2*c2+(c0*c0-THIRD*c2*c1)/c3;
			d=3.0*d/(al*m);
			double t=Math.Abs(d);
			if((t-TOL)<=1.0)
			{
				d=t>1.0?(d>0.0?0.0:Proj.PI):Math.Acos(d);
				lp.phi=Proj.PI*(m*Math.Cos(d*THIRD+PI4_3)-THIRD*c2);
				if(xy.y<0.0) lp.phi=-lp.phi;
				t=r2+TPISQ*(x2-y2+HPISQ);
				lp.lam=Math.Abs(xy.x)<=TOL?0.0:0.5*(r-PISQ+(t<=0.0?0.0:Math.Sqrt(t)))/xy.x;
			}
			else { Proj.pj_ctx_set_errno(ctx, -20); return lp; }

			return lp;
		}

		public override PJ Init()
		{
			es=0.0;
			inv=s_inverse;
			fwd=s_forward;

			return this;
		}
	}
}
