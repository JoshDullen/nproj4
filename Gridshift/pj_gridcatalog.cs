//*****************************************************************************
//
// Project:  PROJ.4
// Purpose:  Code in support of grid catalogs
// Author:   Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2012, Frank Warmerdam <warmerdam@pobox.com>
// Copyright (c) 2014-2015 by the Authors
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

using System.Collections.Generic;
using System.Diagnostics;
using Free.Ports.Proj4.LibcStuff;

namespace Free.Ports.Proj4.Gridshift
{
	public class PJ_GridCatalog
	{
		public struct Entry
		{
			public PJ_Region region;
			public int priority; // higher used before lower
			public double date; // year.fraction
			public string definition; // usually the gridname

			public PJ_GRIDINFO gridinfo;
			public int available; // 0=unknown, 1=true, -1=false
		}

		public string catalog_name;

		public PJ_Region region; // maximum extent of catalog data

		public List<PJ_GridCatalog.Entry> entries;
	}

	public static partial class Grid
	{
		static List<PJ_GridCatalog> grid_catalog_list=new List<PJ_GridCatalog>();

		//**********************************************************************
		//							pj_gc_unloadall()
		//
		//		Deallocate all the grid catalogs (but not the referenced
		//		grids).
		//**********************************************************************
		public static void pj_gc_unloadall(projCtx ctx)
		{
			grid_catalog_list.Clear();
		}

		//**********************************************************************
		//						pj_gc_findcatalog()
		//**********************************************************************
		public static PJ_GridCatalog pj_gc_findcatalog(projCtx ctx, string name)
		{
			lock(gridlock)
			{
				for(int i=0; i<grid_catalog_list.Count; i++)
					if(grid_catalog_list[i].catalog_name==name) return grid_catalog_list[i];
			}

			PJ_GridCatalog catalog=pj_gc_readcatalog(ctx, name);
			if(catalog==null) return null;

			lock(gridlock) grid_catalog_list.Add(catalog);

			return catalog;
		}

		//**********************************************************************
		//						pj_gc_apply_gridshift()
		//**********************************************************************
		public static int pj_gc_apply_gridshift(PJ defn, bool inverse, int point_count, int point_offset, double[] x, double[] y, double[] z)
		{
			int i;

			if(defn.catalog==null)
			{
				defn.catalog=pj_gc_findcatalog(defn.ctx, defn.catalog_name);
				if(defn.catalog==null)
					return defn.ctx.last_errno;
			}

			defn.ctx.last_errno=0;

			for(i=0; i<point_count; i++)
			{
				long io=i*point_offset;

				LP input;
				input.phi=y[io];
				input.lam=x[io];

				// make sure we have appropriate "after" shift file available
				if(defn.last_after_grid==null||
					input.lam<defn.last_after_region.ll_long||input.lam>defn.last_after_region.ur_long||
					input.phi<defn.last_after_region.ll_lat||input.phi>defn.last_after_region.ll_lat)
				{
					defn.last_after_grid=pj_gc_findgrid(defn.ctx, defn.catalog, true, input, defn.datum_date, ref defn.last_after_region, ref defn.last_after_date);
				}

				PJ_GRIDINFO gi=defn.last_after_grid;
				Debug.Assert(gi.child==null);

				// load the grid shift info if we don't have it.
				if(gi.ct.cvs==null&&!pj_gridinfo_load(defn.ctx, gi))
				{
					Proj.pj_ctx_set_errno(defn.ctx, -38);
					return -38;
				}

				LP output_after=nad_cvt(input, inverse, gi.ct);
				if(output_after.lam==Libc.HUGE_VAL)
				{
					if(defn.ctx.debug_level>=PJ_LOG.DEBUG_MAJOR)
					{
						Proj.pj_log(defn.ctx, PJ_LOG.DEBUG_MAJOR, "pj_apply_gridshift(): failed to find a grid shift table for");
						Proj.pj_log(defn.ctx, PJ_LOG.DEBUG_MAJOR, "                      location ({0:F7}dW,{1:F7}dN)", x[io]*Proj.RAD_TO_DEG, y[io]*Proj.RAD_TO_DEG);
					}
					continue;
				}

				if(defn.datum_date==0.0)
				{
					y[io]=output_after.phi;
					x[io]=output_after.lam;
					continue;
				}

				// make sure we have appropriate "before" shift file available
				if(defn.last_before_grid==null||
					input.lam<defn.last_before_region.ll_long||input.lam>defn.last_before_region.ur_long||
					input.phi<defn.last_before_region.ll_lat||input.phi>defn.last_before_region.ll_lat)
				{
					defn.last_before_grid=pj_gc_findgrid(defn.ctx, defn.catalog, false, input, defn.datum_date, ref defn.last_before_region, ref defn.last_before_date);
				}

				gi=defn.last_before_grid;
				Debug.Assert(gi.child==null);

				// load the grid shift info if we don't have it.
				if(gi.ct.cvs==null&&!pj_gridinfo_load(defn.ctx, gi))
				{
					Proj.pj_ctx_set_errno(defn.ctx, -38);
					return -38;
				}

				LP output_before=nad_cvt(input, inverse, gi.ct);
				if(output_before.lam==Libc.HUGE_VAL)
				{
					if(defn.ctx.debug_level>=PJ_LOG.DEBUG_MAJOR)
					{
						Proj.pj_log(defn.ctx, PJ_LOG.DEBUG_MAJOR, "pj_apply_gridshift(): failed to find a grid shift table for");
						Proj.pj_log(defn.ctx, PJ_LOG.DEBUG_MAJOR, "                      location ({0:F7}dW,{1:F7}dN)", x[io]*Proj.RAD_TO_DEG, y[io]*Proj.RAD_TO_DEG);
					}
					continue;
				}

				double mix_ratio=(defn.datum_date-defn.last_before_date)/(defn.last_after_date-defn.last_before_date);

				y[io]=mix_ratio*output_after.phi+(1.0-mix_ratio)*output_before.phi;
				x[io]=mix_ratio*output_after.lam+(1.0-mix_ratio)*output_before.lam;
			}

			return 0;
		}

		//**********************************************************************
		//							pj_c_findgrid()
		//**********************************************************************
		public static PJ_GRIDINFO pj_gc_findgrid(projCtx ctx, PJ_GridCatalog catalog, bool after, LP location, double date, ref PJ_Region optimal_region, ref double grid_date)
		{
			for(int i=0; i<catalog.entries.Count; i++)
			{
				PJ_GridCatalog.Entry entry=catalog.entries[i];

				if((after&&entry.date<date)||(!after&&entry.date>date)) continue;

				if(location.lam<entry.region.ll_long||location.lam>entry.region.ur_long||
					location.phi<entry.region.ll_lat||location.phi>entry.region.ur_lat) continue;

				if(entry.available==-1)
					continue;

				grid_date=entry.date;

				if(entry.gridinfo==null)
				{
					int grid_count;
					PJ_GRIDINFO[] gridlist=pj_gridlist_from_nadgrids(ctx, entry.definition, out grid_count);
					if(grid_count==1) entry.gridinfo=gridlist[0];
				}

				return entry.gridinfo;
			}

			grid_date=0.0;
			optimal_region=new PJ_Region();
			return null;
		}
	}
}
