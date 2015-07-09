//******************************************************************************
// $Id$
//
// Project:	PROJ.4
// Purpose:	Implementation of pj_log() function.
// Author:	Frank Warmerdam, warmerdam@pobox.com
//
//*****************************************************************************
// Copyright (c) 2010, Frank Warmerdam
// Copyright (c) 2008-2011 by the Authors
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

using System;

namespace Free.Ports.Proj4
{
	public static partial class Proj
	{
		//***********************************************************************
		//*							pj_stderr_logger()							*
		//***********************************************************************
		public static void pj_stderr_logger(object app_data, PJ_LOG level, string msg)
		{
			Console.Error.WriteLine(msg);
		}

		//***********************************************************************
		//*								pj_log()								*
		//***********************************************************************
		public static void pj_log(projCtx ctx, PJ_LOG level, string format, params object[] args)
		{
			if(level>ctx.debug_level) return;
			ctx.logger(ctx.app_data, level, string.Format(format, args));
		}
	}
}