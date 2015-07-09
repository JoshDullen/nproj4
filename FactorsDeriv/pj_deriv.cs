using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.FactorsDeriv
{
	public static partial class Factors
	{
		// dervative of (P.fwd) projection
		public static bool pj_deriv(LP lp, double h, PJ P, out DERIVS der)
		{
			der.x_l=der.y_p=der.x_p=der.y_l=0;

			XY t;

			lp.lam+=h;
			lp.phi+=h;
			if(Math.Abs(lp.phi)>Proj.HALFPI) return true;

			h+=h;
			t=P.fwd(lp);
			if(t.x==Libc.HUGE_VAL) return true;

			der.x_l=t.x; der.y_p=t.y; der.x_p=-t.x; der.y_l=-t.y;
			lp.phi-=h;
			if(Math.Abs(lp.phi)>Proj.HALFPI) return true;

			t=P.fwd(lp);
			if(t.x==Libc.HUGE_VAL) return true;

			der.x_l+=t.x; der.y_p-=t.y; der.x_p+=t.x; der.y_l-=t.y;
			lp.lam-=h;

			t=P.fwd(lp);
			if(t.x==Libc.HUGE_VAL) return true;

			der.x_l-=t.x; der.y_p-=t.y; der.x_p+=t.x; der.y_l+=t.y;
			lp.phi+=h;

			t=P.fwd(lp);
			if(t.x==Libc.HUGE_VAL) return true;

			der.x_l-=t.x; der.y_p+=t.y; der.x_p-=t.x; der.y_l+=t.y;
			der.x_l/=(h+=h);
			der.y_p/=h;
			der.x_p/=h;
			der.y_l/=h;

			return false;
		}
	}
}
