#ifndef ISL_INCLUDE_SHADOWCAST_H_
#define ISL_INCLUDE_SHADOWCAST_H_

/* isl_shadowcast - v0.1

   Public domain field of view calcualation using shadowcast technique no
   warranty implied; use at your own risk.

   Do this:
       #define ISL_SHADOWCAST_IMPLEMENTATION
   before you include this file in *one* C or C++ file to create the
   implementation.

   author: Ilya Kolbin (iskolbin@gmail.com)
   url: https://github.com/iskolbin/isl/blob/master/isl_shadowcast.h
   git: git@github.com:iskolbin/isl

   To static link also add:
       #define ISLSC_STATIC

   To adjust stack size (default is 64):
       #define ISLSC_STACK_SIZE 128

   QUICK NOTES:

   This is simple field of view calculation for 2D grid-based games, i.e.
   classic roguelikes. Uses shadowcasting technique with manually converted
   recursion to stack-based implmenentation. Original author of recursive
   shadowcasting is Björn Bergström [bjorn.bergstrom@roguelikedevelopment.org],
   see article on
   [RogueBasin](https://roguebasin.com/index.php/FOV_using_recursive_shadowcasting).

   LICENSE:
     See end of file for license information.
*/

#ifdef ISLSC_STATIC
#define ISLSC_DEF static
#else
#define ISLSC_DEF extern
#endif

#ifndef ISLSC_STACK_SIZE
#define ISLSC_STACK_SIZE 64
#endif

#define ISLSC_OK 0
#define ISLSC_ERROR_STACK_OVERFLOW 1

#ifdef __cplusplus
extern "C" {
#endif

ISLSC_DEF int islsc_update(int x, int y, int radius,
                           void (*update)(int x, int y, void *udata),
                           int (*opaque)(int x, int y, void *udata),
                           void *udata);
#ifdef __cplusplus
}
#endif
#endif // ISL_INCLUDE_SHADOWCAST_H_

#ifdef ISL_SHADOWCAST_IMPLEMENTATION

static int islsc_update(int x, int y, int radius,
                        void (*update)(int x, int y, void *udata),
                        int (*opaque)(int x, int y, void *udata), void *udata) {
  static const int dxys[8][4] = {
      {0, -1, -1, 0}, {-1, 0, 0, -1}, {0, 1, -1, 0}, {1, 0, 0, -1},
      {0, -1, 1, 0},  {-1, 0, 0, 1},  {0, 1, 1, 0},  {1, 0, 0, 1},
  };
  struct islsc_tuple {
    int row;
    float start_slope;
    float end_slope;
  } stack[ISLSC_STACK_SIZE];
  update(x, y, udata);
  for (int i = 0; i < 8; i++) {
    int xx = dxys[i][0], xy = dxys[i][1], yx = dxys[i][2], yy = dxys[i][3];
    stack[0] = (struct islsc_tuple){1, 1.0, 0.0};

    for (int n = 0; n >= 0; n--) {
      int row_ = stack[n].row;
      float start_slope_ = stack[n].start_slope;
      float end_slope_ = stack[n].end_slope;

      if (start_slope_ < end_slope_)
        break;
      float next_start_slope = start_slope_;
      for (int i = row_; i <= radius; i++) {
        int blocked = 0;
        for (int dx = -i, dy = -i; dx <= 0; dx++) {
          float l_slope = (dx - 0.5) / (dy + 0.5);
          float r_slope = (dx + 0.5) / (dy - 0.5);
          if (start_slope_ < r_slope)
            continue;
          else if (end_slope_ > l_slope)
            break;

          int sax = dx * xx + dy * xy;
          int say = dx * yx + dy * yy;
          int abs_sax = sax > 0 ? sax : -sax;
          int abs_say = say > 0 ? say : -say;
          if ((sax < 0 && abs_sax > x) || (say < 0 && abs_say > y))
            continue;

          int ax = x + sax, ay = y + say;
          int radius_sqr = radius * radius;
          if ((dx * dx + dy * dy) < radius_sqr) {
            update(ax, ay, udata);
          }

          int not_visible = opaque(ax, ay, udata);
          if (blocked) {
            if (not_visible) {
              next_start_slope = r_slope;
              continue;
            } else {
              blocked = 0;
              start_slope_ = next_start_slope;
            }
          } else if (not_visible) {
            blocked = 1;
            next_start_slope = r_slope;
            stack[n] = (struct islsc_tuple){i + 1, start_slope_, l_slope};
            n++;
            if (n >= ISLSC_STACK_SIZE)
              return ISLSC_ERROR_STACK_OVERFLOW;
          }
        }
        if (blocked) {
          break;
        }
      }
    }
  }
  return ISLSC_OK;
}

/*
------------------------------------------------------------------------------
This software is available under 2 licenses -- choose whichever you prefer.
------------------------------------------------------------------------------
ALTERNATIVE A - MIT License
Copyright (c) 2025 Ilya Kolbin
Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
------------------------------------------------------------------------------
ALTERNATIVE B - Public Domain (www.unlicense.org)
This is free and unencumbered software released into the public domain.
Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
software, either in source code form or as a compiled binary, for any purpose,
commercial or non-commercial, and by any means.
In jurisdictions that recognize copyright laws, the author or authors of this
software dedicate any and all copyright interest in the software to the public
domain. We make this dedication for the benefit of the public at large and to
the detriment of our heirs and successors. We intend this dedication to be an
overt act of relinquishment in perpetuity of all present and future rights to
this software under copyright law.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
------------------------------------------------------------------------------
*/
#endif // ISL_SHADOWCAST_IMPLEMENTATION
