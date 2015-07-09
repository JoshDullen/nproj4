//******************************************************************************
//
// Project:	PROJ.4
// Purpose:	Implementation of the projCtx thread context object.
// Author:	Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2010, Frank Warmerdam
// Copyright (c) 2008-2015 by the Authors
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
//****************************************************************************

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		static projCtx default_context=new projCtx();
		static volatile bool default_context_initialized=false;

		//***********************************************************************
		//*								pj_get_ctx()							*
		//***********************************************************************
		public static projCtx pj_get_ctx(PJ pj)
		{
			return pj.ctx;
		}

		//***********************************************************************
		//*								pj_set_ctx()							*
		//*																		*
		//*		Note we do not deallocate the old context!						*
		//***********************************************************************
		public static void pj_set_ctx(PJ pj, projCtx ctx)
		{
			pj.ctx=ctx;
		}

		//***********************************************************************
		//*							pj_get_default_ctx()						*
		//***********************************************************************
		public static projCtx pj_get_default_ctx()
		{
			lock(default_context)
			{
				if(!default_context_initialized)
				{
					default_context.last_errno=0;
					default_context.debug_level=PJ_LOG.NONE;
					default_context.logger=pj_stderr_logger;
					default_context.app_data=null;

#if DEBUG

#if PROJ_DEBUG==ERROR
					default_context.debug_level=PJ_LOG.ERROR;
#elif PROJ_DEBUG==DEBUG_MAJOR
					default_context.debug_level=PJ_LOG.DEBUG_MAJOR;
#elif PROJ_DEBUG==DEBUG_MINOR
					default_context.debug_level=PJ_LOG.DEBUG_MINOR;
#else
					default_context.debug_level=PJ_LOG.NONE;
#endif
#endif
					default_context_initialized=true;
				}
			}

			return default_context;
		}

		//***********************************************************************
		//*							pj_ctx_alloc()								*
		//***********************************************************************
		public static projCtx pj_ctx_alloc()
		{
			projCtx ret=new projCtx();
			projCtx ctx=pj_get_default_ctx();
			ret.app_data=ctx.app_data;
			ret.debug_level=ctx.debug_level;
			ret.logger=ctx.logger;

			return ret;
		}

		//***********************************************************************
		//*							pj_ctx_free()								*
		//***********************************************************************
		public static void pj_ctx_free(projCtx ctx)
		{
		}

		//***********************************************************************
		//*							pj_ctx_get_errno()							*
		//***********************************************************************
		public static int pj_ctx_get_errno(projCtx ctx)
		{
			return ctx.last_errno;
		}

		//***********************************************************************
		//*							pj_ctx_set_errno()							*
		//*																		*
		//*		Also sets the global errno.										*
		//***********************************************************************
		public static void pj_ctx_set_errno(projCtx ctx, int errno)
		{
			ctx.last_errno=errno;
			if(errno!=0) pj_errno=errno;
		}

		//***********************************************************************
		//*							pj_ctx_set_debug()							*
		//***********************************************************************
		public static void pj_ctx_set_debug(projCtx ctx, PJ_LOG debug)
		{
			ctx.debug_level=debug;
		}

		//***********************************************************************
		//*							pj_ctx_set_logger()							*
		//***********************************************************************
		public static void pj_ctx_set_logger(projCtx ctx, projCtx.Logger logger)
		{
			ctx.logger=logger;
		}

		//***********************************************************************
		//*							pj_ctx_set_app_data()						*
		//***********************************************************************
		public static void pj_ctx_set_app_data(projCtx ctx, object app_data)
		{
			ctx.app_data=app_data;
		}

		//***********************************************************************
		//*							pj_ctx_get_app_data()						*
		//***********************************************************************
		public static object pj_ctx_get_app_data(projCtx ctx)
		{
			return ctx.app_data;
		}
	}
}
