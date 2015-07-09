using System;

namespace Free.Ports.Proj4.Approximation
{
	public static partial class Approximation
	{
		// sum coefficients less than res
		static void eval(projUV[][] w, int nu, int nv, double res, ref projUV resid)
		{
			resid.u=resid.v=0.0;

			for(int i=0; i<nu; i++)
			{
				for(int j=0; j<nv; j++)
				{
					projUV s=w[i][j];
					double ab=Math.Abs(s.u);
					if(ab<res) resid.u+=ab;
					ab=Math.Abs(s.v);
					if(ab<res) resid.v+=ab;
				}
			}
		}

		// create power series structure
		static Tseries makeT(int nru, int nrv)
		{
			try
			{
				Tseries T=new Tseries();
				T.cu=new PW_COEF[nru];
				T.cv=new PW_COEF[nrv];
				for(int i=0; i<nru; i++) T.cu[i].c=null;
				for(int i=0; i<nrv; i++) T.cv[i].c=null;
				return T;
			}
			catch
			{
				return null;
			}
		}

		public static Tseries mk_cheby(projUV a, projUV b, double res, ref projUV resid, ProjFunc func, int nu, int nv, bool power)
		{
			try
			{
				projUV[][] w=new projUV[nu][];
				for(int k=0; k<nu; k++) w[k]=new projUV[nv];

				if(bchgen(a, b, nu, nv, w, func)) return null;

				// analyse coefficients and adjust until residual OK
				double cutres=res;
				int i=4;
				for(; i>0; i--)
				{
					eval(w, nu, nv, cutres, ref resid);
					if(resid.u<res&&resid.v<res) break;
					cutres*=0.5;
				}

				// warn of too many tries 
				if(i<=0) resid.u=-resid.u;

				double ab;
				double[] p;

				int[] ncu=new int[nu];
				int[] ncv=new int[nv];

				// apply cut resolution and set pointers
				int nru=0, nrv=0;
				for(int j=0; j<nu; ++j)
				{
					ncu[j]=ncv[j]=0; // clear column maxes

					projUV[] s=w[j];
					for(i=0; i<nv; i++)
					{
						// < resolution ?
						ab=Math.Abs(s[i].u);
						if(ab<cutres) s[i].u=0.0;	// clear coefficient
						else ncu[j]=i+1;			// update column max

						ab=Math.Abs(s[i].v);
						if(ab<cutres) s[i].v=0.0;	// same for v coef's
						else ncv[j]=i+1;
					}
					if(ncu[j]!=0) nru=j+1;			// update row max
					if(ncv[j]!=0) nrv=j+1;
				}

				if(power)
				{ // convert to bivariate power series
					if(!bch2bps(a, b, w, nu, nv)) return null;

					// possible change in some row counts, so readjust
					nru=nrv=0;
					for(int j=0; j<nu; ++j)
					{
						ncu[j]=ncv[j]=0; // clear column maxes

						projUV[] s=w[j];

						for(i=0; i<nv; i++)
						{
							if(s[i].u!=0) ncu[j]=i+1;	// update column max
							if(s[i].v!=0) ncv[j]=i+1;
						}
						if(ncu[j]!=0) nru=j+1;			// update row max
						if(ncv[j]!=0) nrv=j+1;
					}

					Tseries T=makeT(nru, nrv);
					if(T!=null)
					{
						T.a=a;
						T.b=b;
						T.mu=nru-1;
						T.mv=nrv-1;
						T.power=1;

						for(i=0; i<nru; ++i) // store coefficient rows for u
						{
							T.cu[i].m=ncu[i];
							if(T.cu[i].m!=0)
							{
								p=T.cu[i].c=new double[ncu[i]];
								projUV[] s=w[i];
								for(int j=0; j<ncu[i]; j++) p[j]=s[j].u;
							}
						}

						for(i=0; i<nrv; ++i) // same for v
						{
							T.cv[i].m=ncv[i];
							if(T.cv[i].m!=0)
							{
								p=T.cv[i].c=new double[ncv[i]];
								projUV[] s=w[i];
								for(int j=0; j<ncv[i]; j++) p[j]=s[j].v;
							}
						}
					}
					return T;
				}
				else
				{
					Tseries T=makeT(nru, nrv);
					if(T!=null)
					{
						// else make returned Chebyshev coefficient structure
						T.mu=nru-1; // save row degree
						T.mv=nrv-1;
						T.a.u=a.u+b.u; // set argument scaling
						T.a.v=a.v+b.v;
						T.b.u=1.0/(b.u-a.u);
						T.b.v=1.0/(b.v-a.v);
						T.power=0;

						for(i=0; i<nru; ++i) // store coefficient rows for u
						{
							T.cu[i].m=ncu[i];
							if(T.cu[i].m!=0)
							{
								p=T.cu[i].c=new double[ncu[i]];
								projUV[] s=w[i];
								for(int j=0; j<ncu[i]; j++) p[j]=s[j].u;
							}
						}

						for(i=0; i<nrv; ++i) // same for v
						{
							T.cv[i].m=ncv[i];
							if(T.cv[i].m!=0)
							{
								p=T.cv[i].c=new double[ncv[i]];
								projUV[] s=w[i];
								for(int j=0; j<ncv[i]; j++) p[j]=s[j].v;
							}
						}
					}
					return T;
				}
			}
			catch
			{
				return null;
			}
		}
	}
}
