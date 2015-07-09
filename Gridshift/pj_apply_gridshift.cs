//*****************************************************************************
//
// Project:	PROJ.4
// Purpose:	Apply datum shifts based on grid shift files (normally NAD27 to
//			NAD83 or the reverse). This module is responsible for keeping
//			a list of loaded grids, and calling with each one that is 
//			allowed for a given datum (expressed as the nadgrids= parameter).
// Author:	Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2000, Frank Warmerdam <warmerdam@pobox.com>
// Copyright (c) 2009-2015 by the Authors
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//*****************************************************************************

using System;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Gridshift
{
	public static partial class Grid
	{
		//**********************************************************************
		//							pj_apply_gridshift()
		//
		//		This is the externally callable interface - part of the
		//		public API - though it is not used internally any more and I
		//		doubt it is used by any other applications. But we preserve
		//		it to honour our public api.
		//**********************************************************************
		public static int pj_apply_gridshift(projCtx ctx, string nadgrids, bool inverse, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			int grid_count=0;
			PJ_GRIDINFO[] gridlist=pj_gridlist_from_nadgrids(ctx, nadgrids, out grid_count);
			if(gridlist==null||grid_count==0) return ctx.last_errno;

			return pj_apply_gridshift_3(ctx, gridlist, grid_count, inverse, point_count, point_offset, x, y, z);
		}

		//***********************************************************************
		//							pj_apply_gridshift_2()
		//
		//		This implmentation takes uses the gridlist from a coordinate
		//		system definition.  If the gridlist has not yet been
		//		populated in the coordinate system definition we set it up
		//		now.
		//***********************************************************************
		public static int pj_apply_gridshift_2(PJ defn, bool inverse, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			if(defn.catalog_name!=null) return pj_gc_apply_gridshift(defn, inverse, point_count, point_offset, x, y, z);

			if(defn.gridlist==null)
			{
				defn.gridlist=pj_gridlist_from_nadgrids(Proj.pj_get_ctx(defn), Proj.pj_param_s(defn.ctx, defn.parameters, "nadgrids"), out defn.gridlist_count);

				if(defn.gridlist==null||defn.gridlist_count==0)
					return defn.ctx.last_errno;
			}

			return pj_apply_gridshift_3(Proj.pj_get_ctx(defn), defn.gridlist, defn.gridlist_count, inverse, point_count, point_offset, x, y, z);
		}

		//**********************************************************************
		//						pj_apply_gridshift_3()
		//
		//		This is the real workhorse, given a gridlist.
		//**********************************************************************
		static int debug_count_apply_gridshift_3=0; // TODO
		public static int pj_apply_gridshift_3(projCtx ctx, PJ_GRIDINFO[] tables, int grid_count, bool inverse, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			if(tables==null||grid_count==0)
			{
				Proj.pj_ctx_set_errno(ctx, -38);
				return -38;
			}

			ctx.last_errno=0;

			for(int i=0; i<point_count; i++)
			{
				long io=i*point_offset;
				LP input, output;
				int itable;

				input.phi=y[io];
				input.lam=x[io];
				output.phi=Libc.HUGE_VAL;
				output.lam=Libc.HUGE_VAL;

				// keep trying till we find a table that works
				for(itable=0; itable<grid_count; itable++)
				{
					PJ_GRIDINFO gi=tables[itable];
					CTABLE ct=gi.ct;

					double epsilon=(Math.Abs(ct.del.phi)+Math.Abs(ct.del.lam))/10000.0;

					// skip tables that don't match our point at all.
					if(ct.ll.phi-epsilon>input.phi||
						ct.ll.lam-epsilon>input.lam||
						ct.ll.phi+(ct.lim.phi-1)*ct.del.phi+epsilon<input.phi||
						ct.ll.lam+(ct.lim.lam-1)*ct.del.lam+epsilon<input.lam) continue;

					// If we have child nodes, check to see if any of them apply.
					while(gi.child!=null)
					{
						PJ_GRIDINFO child=gi.child;

						for(; child!=null; child=child.next)
						{
							CTABLE ct1=child.ct;
							if(ct1.ll.phi-epsilon>input.phi||
								ct1.ll.lam-epsilon>input.lam||
								ct1.ll.phi+(ct1.lim.phi-1)*ct1.del.phi+epsilon<input.phi||
								ct1.ll.lam+(ct1.lim.lam-1)*ct1.del.lam+epsilon<input.lam) continue;

							break;
						}

						// If we didn't find a child then nothing more to do,
						if(child==null) break;

						// otherwise use the child, first checking it's children.
						gi=child;
						ct=child.ct;
					}

					// load the grid shift info if we don't have it.
					if(ct.cvs==null&&!pj_gridinfo_load(ctx, gi))
					{
						Proj.pj_ctx_set_errno(ctx, -38);
						return -38;
					}

					output=nad_cvt(input, inverse, ct);
					if(output.lam!=Libc.HUGE_VAL)
					{
						if(debug_count_apply_gridshift_3++<20) Proj.pj_log(ctx, PJ_LOG.DEBUG_MINOR, "pj_apply_gridshift(): used {0}", ct.id);
						break;
					}
				}

				if(output.lam==Libc.HUGE_VAL)
				{
					if(ctx.debug_level>=PJ_LOG.DEBUG_MAJOR)
					{
						Proj.pj_log(ctx, PJ_LOG.DEBUG_MAJOR, "pj_apply_gridshift(): failed to find a grid shift table for\n\t\t\tlocation ({0}dW,{1}dN)", x[io]*Proj.RAD_TO_DEG, y[io]*Proj.RAD_TO_DEG);
						for(itable=0; itable<grid_count; itable++)
						{
							PJ_GRIDINFO gi=tables[itable];
							if(itable==0) Proj.pj_log(ctx, PJ_LOG.DEBUG_MAJOR, "   tried: {0}", gi.gridname);
							else Proj.pj_log(ctx, PJ_LOG.DEBUG_MAJOR, ",{0}", gi.gridname);
						}
					}

					// We don't actually have any machinery currently to set the
					// following macro, so this is mostly kept here to make it clear
					// how we ought to operate if we wanted to make it super clear
					// that an error has occured when points are outside our available
					// datum shift areas. But if this is on, we will find that "low
					// value" points on the fringes of some datasets will completely
					// fail causing lots of problems when it is more or less ok to
					// just not apply a datum shift. So rather than deal with
					// that we just fallback to no shift. (see also bug #45).
#if ERR_GRID_AREA_TRANSIENT_SEVERE
					y[io]=Libc.HUGE_VAL;
					x[io]=Libc.HUGE_VAL;
#else
					// leave x/y unshifted.
#endif
				}
				else
				{
					y[io]=output.phi;
					x[io]=output.lam;
				}
			}

			return 0;
		}
	}
}
