using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Gridshift
{
	public static partial class Grid
	{
		const double TOL12=1.0e-12;

		public static LP nad_cvt(LP @in, bool inverse, CTABLE ct)
		{
			const int MAX_TRY=9;

			if(@in.lam==Libc.HUGE_VAL) return @in;

			// normalize input to ll origin
			LP tb=@in;
			tb.lam-=ct.ll.lam;
			tb.phi-=ct.ll.phi;
			tb.lam=Proj.adjlon(tb.lam-Proj.PI)+Proj.PI;
			LP t=nad_intr(tb, ct);
			if(inverse)
			{
				LP del, dif;
				int i=MAX_TRY;

				if(t.lam==Libc.HUGE_VAL) return t;
				t.lam=tb.lam+t.lam;
				t.phi=tb.phi-t.phi;

				do
				{
					del=nad_intr(t, ct);

					// This case used to return failure, but I have
					// changed it to return the first order approximation
					// of the inverse shift. This avoids cases where the
					// grid shift *into* this grid came from another grid.
					// While we aren't returning optimally correct results
					// I feel a close result in this case is better than
					// no result. NFW
					// To demonstrate use -112.5839956 49.4914451 against
					// the NTv2 grid shift file from Canada.
					if(del.lam==Libc.HUGE_VAL)
					{
#if DEBUG
						Console.Error.WriteLine("Inverse grid shift iteration failed, presumably at grid edge.\nUsing first approximation.");
#endif
						break;
					}

					t.lam-=dif.lam=t.lam-del.lam-tb.lam;
					t.phi-=dif.phi=t.phi+del.phi-tb.phi;
				} while((i--)!=0&&Math.Abs(dif.lam)>TOL12&&Math.Abs(dif.phi)>TOL12);

				if(i<0)
				{
#if DEBUG
					Console.Error.WriteLine("Inverse grid shift iterator failed to converge.");
#endif
					t.lam=t.phi=Libc.HUGE_VAL;
					return t;
				}

				@in.lam=Proj.adjlon(t.lam+ct.ll.lam);
				@in.phi=t.phi+ct.ll.phi;
			}
			else
			{
				if(t.lam==Libc.HUGE_VAL) @in=t;
				else
				{
					@in.lam-=t.lam;
					@in.phi+=t.phi;
				}
			}
			return @in;
		}
	}
}
