int domain_init_turb(void)
{
  int i, j, k;    // iterator
  int C, W, E, S, N, B, T;
  real tmp;

  // make sure there are enough GPU devices in the given range
  if(nsubdom > dev_end - dev_start + 1) {
    return EXIT_FAILURE;
  }

  // calculate domain sizes
  Dom.xl = Dom.xe - Dom.xs;
  Dom.yl = Dom.ye - Dom.ys;
  Dom.zl = Dom.ze - Dom.zs;

  // calculate cell sizes
  Dom.dx = Dom.xl / Dom.xn;
  Dom.dy = Dom.yl / Dom.yn;
  Dom.dz = Dom.zl / Dom.zn;

  // set up grids
  // Gcc
  Dom.Gcc.is = DOM_BUF;
  Dom.Gcc.isb = Dom.Gcc.is - DOM_BUF;
  Dom.Gcc.in = Dom.xn;
  Dom.Gcc.inb = Dom.Gcc.in + 2 * DOM_BUF;
  Dom.Gcc.ie = Dom.Gcc.is + Dom.Gcc.in;
  Dom.Gcc.ieb = Dom.Gcc.ie + DOM_BUF;

  Dom.Gcc.js = DOM_BUF;
  Dom.Gcc.jsb = Dom.Gcc.js - DOM_BUF;
  Dom.Gcc.jn = Dom.yn;
  Dom.Gcc.jnb = Dom.Gcc.jn + 2 * DOM_BUF;
  Dom.Gcc.je = Dom.Gcc.js + Dom.Gcc.jn;
  Dom.Gcc.jeb = Dom.Gcc.je + DOM_BUF;

  Dom.Gcc.ks = DOM_BUF;
  Dom.Gcc.ksb = Dom.Gcc.ks - DOM_BUF;
  Dom.Gcc.kn = Dom.zn;
  Dom.Gcc.knb = Dom.Gcc.kn + 2 * DOM_BUF;
  Dom.Gcc.ke = DOM_BUF + Dom.Gcc.kn;
  Dom.Gcc.keb = Dom.Gcc.ke + DOM_BUF;

  Dom.Gcc.s1 = Dom.Gcc.in;
  Dom.Gcc.s2 = Dom.Gcc.s1 * Dom.Gcc.jn;
  Dom.Gcc.s3 = Dom.Gcc.s2 * Dom.Gcc.kn;
  Dom.Gcc.s1b = Dom.Gcc.inb;
  Dom.Gcc.s2b = Dom.Gcc.s1b * Dom.Gcc.jnb;
  Dom.Gcc.s3b = Dom.Gcc.s2b * Dom.Gcc.knb;

  // Gfx
  Dom.Gfx.is = DOM_BUF;
  Dom.Gfx.isb = Dom.Gfx.is - DOM_BUF;
  Dom.Gfx.in = Dom.xn + 1;
  Dom.Gfx.inb = Dom.Gfx.in + 2 * DOM_BUF;
  Dom.Gfx.ie = Dom.Gfx.is + Dom.Gfx.in;
  Dom.Gfx.ieb = Dom.Gfx.ie + DOM_BUF;

  Dom.Gfx.js = DOM_BUF;
  Dom.Gfx.jsb = Dom.Gfx.js - DOM_BUF;
  Dom.Gfx.jn = Dom.yn;
  Dom.Gfx.jnb = Dom.Gfx.jn + 2 * DOM_BUF;
  Dom.Gfx.je = Dom.Gfx.js + Dom.Gfx.jn;
  Dom.Gfx.jeb = Dom.Gfx.je + DOM_BUF;

  Dom.Gfx.ks = DOM_BUF;
  Dom.Gfx.ksb = Dom.Gfx.ks - DOM_BUF;
  Dom.Gfx.kn = Dom.zn;
  Dom.Gfx.knb = Dom.Gfx.kn + 2 * DOM_BUF;
  Dom.Gfx.ke = Dom.Gfx.ks + Dom.Gfx.kn;
  Dom.Gfx.keb = Dom.Gfx.ke + DOM_BUF;

  Dom.Gfx.s1 = Dom.Gfx.in;
  Dom.Gfx.s2 = Dom.Gfx.s1 * Dom.Gfx.jn;
  Dom.Gfx.s3 = Dom.Gfx.s2 * Dom.Gfx.kn;
  Dom.Gfx.s1b = Dom.Gfx.inb;
  Dom.Gfx.s2b = Dom.Gfx.s1b * Dom.Gfx.jnb;
  Dom.Gfx.s3b = Dom.Gfx.s2b * Dom.Gfx.knb;

  // Gfy
  Dom.Gfy.is = DOM_BUF;
  Dom.Gfy.isb = Dom.Gfy.is - DOM_BUF;
  Dom.Gfy.in = Dom.xn;
  Dom.Gfy.inb = Dom.Gfy.in + 2 * DOM_BUF;
  Dom.Gfy.ie = Dom.Gfy.is + Dom.Gfy.in;
  Dom.Gfy.ieb = Dom.Gfy.ie + DOM_BUF;

  Dom.Gfy.js = DOM_BUF;
  Dom.Gfy.jsb = Dom.Gfy.js - DOM_BUF;
  Dom.Gfy.jn = Dom.yn + 1;
  Dom.Gfy.jnb = Dom.Gfy.jn + 2 * DOM_BUF;
  Dom.Gfy.je = Dom.Gfy.js + Dom.Gfy.jn;
  Dom.Gfy.jeb = Dom.Gfy.je + DOM_BUF;

  Dom.Gfy.ks = DOM_BUF;
  Dom.Gfy.ksb = Dom.Gfy.ks - DOM_BUF;
  Dom.Gfy.kn = Dom.zn;
  Dom.Gfy.knb = Dom.Gfy.kn + 2 * DOM_BUF;
  Dom.Gfy.ke = Dom.Gfy.ks + Dom.Gfy.kn;
  Dom.Gfy.keb = Dom.Gfy.ke + DOM_BUF;

  Dom.Gfy.s1 = Dom.Gfy.in;
  Dom.Gfy.s2 = Dom.Gfy.s1 * Dom.Gfy.jn;
  Dom.Gfy.s3 = Dom.Gfy.s2 * Dom.Gfy.kn;
  Dom.Gfy.s1b = Dom.Gfy.inb;
  Dom.Gfy.s2b = Dom.Gfy.s1b * Dom.Gfy.jnb;
  Dom.Gfy.s3b = Dom.Gfy.s2b * Dom.Gfy.knb;

  // Gfz
  Dom.Gfz.is = DOM_BUF;
  Dom.Gfz.isb = Dom.Gfz.is - DOM_BUF;
  Dom.Gfz.in = Dom.xn;
  Dom.Gfz.inb = Dom.Gfz.in + 2 * DOM_BUF;
  Dom.Gfz.ie = Dom.Gfz.is + Dom.Gfz.in;
  Dom.Gfz.ieb = Dom.Gfz.ie + DOM_BUF;

  Dom.Gfz.js = DOM_BUF;
  Dom.Gfz.jsb = Dom.Gfz.js - DOM_BUF;
  Dom.Gfz.jn = Dom.yn;
  Dom.Gfz.jnb = Dom.Gfz.jn + 2 * DOM_BUF;
  Dom.Gfz.je = Dom.Gfz.js + Dom.Gfz.jn;
  Dom.Gfz.jeb = Dom.Gfz.je + DOM_BUF;

  Dom.Gfz.ks = DOM_BUF;
  Dom.Gfz.ksb = Dom.Gfz.ks - DOM_BUF;
  Dom.Gfz.kn = Dom.zn + 1;
  Dom.Gfz.knb = Dom.Gfz.kn + 2 * DOM_BUF;
  Dom.Gfz.ke = Dom.Gfz.ks + Dom.Gfz.kn;
  Dom.Gfz.keb = Dom.Gfz.ke + DOM_BUF;

  Dom.Gfz.s1 = Dom.Gfz.in;
  Dom.Gfz.s2 = Dom.Gfz.s1 * Dom.Gfz.jn;
  Dom.Gfz.s3 = Dom.Gfz.s2 * Dom.Gfz.kn;
  Dom.Gfz.s1b = Dom.Gfz.inb;
  Dom.Gfz.s2b = Dom.Gfz.s1b * Dom.Gfz.jnb;
  Dom.Gfz.s3b = Dom.Gfz.s2b * Dom.Gfz.knb;

  // initialize subdomains
  for(i = 0; i < nsubdom; i++) {
    dom[i].xl = dom[i].xe - dom[i].xs;
    dom[i].yl = dom[i].ye - dom[i].ys;
    dom[i].zl = dom[i].ze - dom[i].zs;
    dom[i].dx = dom[i].xl / dom[i].xn;
    dom[i].dy = dom[i].yl / dom[i].yn;
    dom[i].dz = dom[i].zl / dom[i].zn;

    // TODO: this algorithm will fail if subdomains are not numbered in
    // increasing order.  This will need to be fixed before going to a machine
    // with more than two GPU's.
    // Gcc
    if(dom[i].W > -1)
      dom[i].Gcc.is = dom[dom[i].W].Gcc.ie;
    else
      dom[i].Gcc.is = DOM_BUF;
    dom[i].Gcc.isb = dom[i].Gcc.is - DOM_BUF;
    dom[i].Gcc.in = dom[i].xn;
    dom[i].Gcc.inb = dom[i].Gcc.in + 2 * DOM_BUF;
    dom[i].Gcc.ie = dom[i].Gcc.is + dom[i].Gcc.in;
    dom[i].Gcc.ieb = dom[i].Gcc.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gcc.js = dom[dom[i].S].Gcc.je;
    else
      dom[i].Gcc.js = DOM_BUF;
    dom[i].Gcc.jsb = dom[i].Gcc.js - DOM_BUF;
    dom[i].Gcc.jn = dom[i].yn;
    dom[i].Gcc.jnb = dom[i].Gcc.jn + 2 * DOM_BUF;
    dom[i].Gcc.je = dom[i].Gcc.js + dom[i].Gcc.jn;
    dom[i].Gcc.jeb = dom[i].Gcc.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gcc.ks = dom[dom[i].B].Gcc.ke;
    else
      dom[i].Gcc.ks = DOM_BUF;
    dom[i].Gcc.ksb = dom[i].Gcc.ks - DOM_BUF;
    dom[i].Gcc.kn = dom[i].zn;
    dom[i].Gcc.knb = dom[i].Gcc.kn + 2 * DOM_BUF;
    dom[i].Gcc.ke = DOM_BUF + dom[i].Gcc.kn;
    dom[i].Gcc.keb = dom[i].Gcc.ke + DOM_BUF;

    dom[i].Gcc.s1 = dom[i].Gcc.in;
    dom[i].Gcc.s2 = dom[i].Gcc.s1 * dom[i].Gcc.jn;
    dom[i].Gcc.s3 = dom[i].Gcc.s2 * dom[i].Gcc.kn;
    dom[i].Gcc.s1b = dom[i].Gcc.inb;
    dom[i].Gcc.s2b = dom[i].Gcc.s1b * dom[i].Gcc.jnb;
    dom[i].Gcc.s3b = dom[i].Gcc.s2b * dom[i].Gcc.knb;

    dom[i].Gcc._is = DOM_BUF;
    dom[i].Gcc._isb = dom[i].Gcc._is - DOM_BUF;
    dom[i].Gcc._in = dom[i].xn;
    dom[i].Gcc._inb = dom[i].Gcc._in + 2 * DOM_BUF;
    dom[i].Gcc._ie = dom[i].Gcc._is + dom[i].Gcc._in;
    dom[i].Gcc._ieb = dom[i].Gcc._ie + DOM_BUF;

    dom[i].Gcc._js = DOM_BUF;
    dom[i].Gcc._jsb = dom[i].Gcc._js - DOM_BUF;
    dom[i].Gcc._jn = dom[i].yn;
    dom[i].Gcc._jnb = dom[i].Gcc._jn + 2 * DOM_BUF;
    dom[i].Gcc._je = dom[i].Gcc._js + dom[i].Gcc._jn;
    dom[i].Gcc._jeb = dom[i].Gcc._je + DOM_BUF;

    dom[i].Gcc._ks = DOM_BUF;
    dom[i].Gcc._ksb = dom[i].Gcc._ks - DOM_BUF;
    dom[i].Gcc._kn = dom[i].zn;
    dom[i].Gcc._knb = dom[i].Gcc._kn + 2 * DOM_BUF;
    dom[i].Gcc._ke = dom[i].Gcc._ks + dom[i].Gcc._kn;
    dom[i].Gcc._keb = dom[i].Gcc._ke + DOM_BUF;

    dom[i].Gcc._s1 = dom[i].Gcc._in;
    dom[i].Gcc._s2 = dom[i].Gcc._s1 * dom[i].Gcc._jn;
    dom[i].Gcc._s3 = dom[i].Gcc._s2 * dom[i].Gcc._kn;
    dom[i].Gcc._s1b = dom[i].Gcc._inb;
    dom[i].Gcc._s2b = dom[i].Gcc._s1b * dom[i].Gcc._jnb;
    dom[i].Gcc._s3b = dom[i].Gcc._s2b * dom[i].Gcc._knb;

    // Gfx
    if(dom[i].W > -1)
      dom[i].Gfx.is = dom[dom[i].W].Gfx.ie - 1;
    else
      dom[i].Gfx.is = DOM_BUF;
    dom[i].Gfx.isb = dom[i].Gfx.is - DOM_BUF;
    dom[i].Gfx.in = dom[i].xn + 1;
    dom[i].Gfx.inb = dom[i].Gfx.in + 2 * DOM_BUF;
    dom[i].Gfx.ie = dom[i].Gfx.is + dom[i].Gfx.in;
    dom[i].Gfx.ieb = dom[i].Gfx.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfx.js = dom[dom[i].S].Gfx.je - 1;
    else
      dom[i].Gfx.js = DOM_BUF;
    dom[i].Gfx.jsb = dom[i].Gfx.js - DOM_BUF;
    dom[i].Gfx.jn = dom[i].yn;
    dom[i].Gfx.jnb = dom[i].Gfx.jn + 2 * DOM_BUF;
    dom[i].Gfx.je = dom[i].Gfx.js + dom[i].Gfx.jn;
    dom[i].Gfx.jeb = dom[i].Gfx.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfx.ks = dom[dom[i].B].Gfx.ke - 1;
    else
      dom[i].Gfx.ks = DOM_BUF;
    dom[i].Gfx.ksb = dom[i].Gfx.ks - DOM_BUF;
    dom[i].Gfx.kn = dom[i].zn;
    dom[i].Gfx.knb = dom[i].Gfx.kn + 2 * DOM_BUF;
    dom[i].Gfx.ke = dom[i].Gfx.ks + dom[i].Gfx.kn;
    dom[i].Gfx.keb = dom[i].Gfx.ke + DOM_BUF;

    dom[i].Gfx.s1 = dom[i].Gfx.in;
    dom[i].Gfx.s2 = dom[i].Gfx.s1 * dom[i].Gfx.jn;
    dom[i].Gfx.s3 = dom[i].Gfx.s2 * dom[i].Gfx.kn;
    dom[i].Gfx.s1b = dom[i].Gfx.inb;
    dom[i].Gfx.s2b = dom[i].Gfx.s1b * dom[i].Gfx.jnb;
    dom[i].Gfx.s3b = dom[i].Gfx.s2b * dom[i].Gfx.knb;

    dom[i].Gfx._is = DOM_BUF;
    dom[i].Gfx._isb = dom[i].Gfx._is - DOM_BUF;
    dom[i].Gfx._in = dom[i].xn + 1;
    dom[i].Gfx._inb = dom[i].Gfx._in + 2 * DOM_BUF;
    dom[i].Gfx._ie = dom[i].Gfx._is + dom[i].Gfx._in;
    dom[i].Gfx._ieb = dom[i].Gfx._ie + DOM_BUF;

    dom[i].Gfx._js = DOM_BUF;
    dom[i].Gfx._jsb = dom[i].Gfx._js - DOM_BUF;
    dom[i].Gfx._jn = dom[i].yn;
    dom[i].Gfx._jnb = dom[i].Gfx._jn + 2 * DOM_BUF;
    dom[i].Gfx._je = dom[i].Gfx._js + dom[i].Gfx._jn;
    dom[i].Gfx._jeb = dom[i].Gfx._je + DOM_BUF;

    dom[i].Gfx._ks = DOM_BUF;
    dom[i].Gfx._ksb = dom[i].Gfx._ks - DOM_BUF;
    dom[i].Gfx._kn = dom[i].zn;
    dom[i].Gfx._knb = dom[i].Gfx._kn + 2 * DOM_BUF;
    dom[i].Gfx._ke = dom[i].Gfx._ks + dom[i].Gfx._kn;
    dom[i].Gfx._keb = dom[i].Gfx._ke + DOM_BUF;

    dom[i].Gfx._s1 = dom[i].Gfx._in;
    dom[i].Gfx._s2 = dom[i].Gfx._s1 * dom[i].Gfx._jn;
    dom[i].Gfx._s3 = dom[i].Gfx._s2 * dom[i].Gfx._kn;
    dom[i].Gfx._s1b = dom[i].Gfx._inb;
    dom[i].Gfx._s2b = dom[i].Gfx._s1b * dom[i].Gfx._jnb;
    dom[i].Gfx._s3b = dom[i].Gfx._s2b * dom[i].Gfx._knb;

    // Gfy
    if(dom[i].W > -1)
      dom[i].Gfy.is = dom[dom[i].W].Gfy.ie;
    else
      dom[i].Gfy.is = DOM_BUF;
    dom[i].Gfy.isb = dom[i].Gfy.is - DOM_BUF;
    dom[i].Gfy.in = dom[i].xn;
    dom[i].Gfy.inb = dom[i].Gfy.in + 2 * DOM_BUF;
    dom[i].Gfy.ie = dom[i].Gfy.is + dom[i].Gfy.in;
    dom[i].Gfy.ieb = dom[i].Gfy.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfy.js = dom[dom[i].S].Gfy.je;
    else
      dom[i].Gfy.js = DOM_BUF;
    dom[i].Gfy.jsb = dom[i].Gfy.js - DOM_BUF;
    dom[i].Gfy.jn = dom[i].yn + 1;
    dom[i].Gfy.jnb = dom[i].Gfy.jn + 2 * DOM_BUF;
    dom[i].Gfy.je = dom[i].Gfy.js + dom[i].Gfy.jn;
    dom[i].Gfy.jeb = dom[i].Gfy.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfy.ks = dom[dom[i].B].Gfy.ke;
    else
      dom[i].Gfy.ks = DOM_BUF;
    dom[i].Gfy.ksb = dom[i].Gfy.ks - DOM_BUF;
    dom[i].Gfy.kn = dom[i].zn;
    dom[i].Gfy.knb = dom[i].Gfy.kn + 2 * DOM_BUF;
    dom[i].Gfy.ke = dom[i].Gfy.ks + dom[i].Gfy.kn;
    dom[i].Gfy.keb = dom[i].Gfy.ke + DOM_BUF;

    dom[i].Gfy.s1 = dom[i].Gfy.in;
    dom[i].Gfy.s2 = dom[i].Gfy.s1 * dom[i].Gfy.jn;
    dom[i].Gfy.s3 = dom[i].Gfy.s2 * dom[i].Gfy.kn;
    dom[i].Gfy.s1b = dom[i].Gfy.inb;
    dom[i].Gfy.s2b = dom[i].Gfy.s1b * dom[i].Gfy.jnb;
    dom[i].Gfy.s3b = dom[i].Gfy.s2b * dom[i].Gfy.knb;

    dom[i].Gfy._is = DOM_BUF;
    dom[i].Gfy._isb = dom[i].Gfy._is - DOM_BUF;
    dom[i].Gfy._in = dom[i].xn;
    dom[i].Gfy._inb = dom[i].Gfy._in + 2 * DOM_BUF;
    dom[i].Gfy._ie = dom[i].Gfy._is + dom[i].Gfy._in;
    dom[i].Gfy._ieb = dom[i].Gfy._ie + DOM_BUF;

    dom[i].Gfy._js = DOM_BUF;
    dom[i].Gfy._jsb = dom[i].Gfy._js - DOM_BUF;
    dom[i].Gfy._jn = dom[i].yn + 1;
    dom[i].Gfy._jnb = dom[i].Gfy._jn + 2 * DOM_BUF;
    dom[i].Gfy._je = dom[i].Gfy._js + dom[i].Gfy._jn;
    dom[i].Gfy._jeb = dom[i].Gfy._je + DOM_BUF;

    dom[i].Gfy._ks = DOM_BUF;
    dom[i].Gfy._ksb = dom[i].Gfy._ks - DOM_BUF;
    dom[i].Gfy._kn = dom[i].zn;
    dom[i].Gfy._knb = dom[i].Gfy._kn + 2 * DOM_BUF;
    dom[i].Gfy._ke = dom[i].Gfy._ks + dom[i].Gfy._kn;
    dom[i].Gfy._keb = dom[i].Gfy._ke + DOM_BUF;

    dom[i].Gfy._s1 = dom[i].Gfy._in;
    dom[i].Gfy._s2 = dom[i].Gfy._s1 * dom[i].Gfy._jn;
    dom[i].Gfy._s3 = dom[i].Gfy._s2 * dom[i].Gfy._kn;
    dom[i].Gfy._s1b = dom[i].Gfy._inb;
    dom[i].Gfy._s2b = dom[i].Gfy._s1b * dom[i].Gfy._jnb;
    dom[i].Gfy._s3b = dom[i].Gfy._s2b * dom[i].Gfy._knb;

    // Gfz
    if(dom[i].W > -1)
      dom[i].Gfz.is = dom[dom[i].W].Gfz.ie;
    else
      dom[i].Gfz.is = DOM_BUF;
    dom[i].Gfz.isb = dom[i].Gfz.is - DOM_BUF;
    dom[i].Gfz.in = dom[i].xn;
    dom[i].Gfz.inb = dom[i].Gfz.in + 2 * DOM_BUF;
    dom[i].Gfz.ie = dom[i].Gfz.is + dom[i].Gfz.in;
    dom[i].Gfz.ieb = dom[i].Gfz.ie + DOM_BUF;

    if(dom[i].S > -1)
      dom[i].Gfz.js = dom[dom[i].S].Gfz.je;
    else
      dom[i].Gfz.js = DOM_BUF;
    dom[i].Gfz.jsb = dom[i].Gfz.js - DOM_BUF;
    dom[i].Gfz.jn = dom[i].yn;
    dom[i].Gfz.jnb = dom[i].Gfz.jn + 2 * DOM_BUF;
    dom[i].Gfz.je = dom[i].Gfz.js + dom[i].Gfz.jn;
    dom[i].Gfz.jeb = dom[i].Gfz.je + DOM_BUF;

    if(dom[i].B > -1)
      dom[i].Gfz.ks = dom[dom[i].B].Gfz.ke;
    else
      dom[i].Gfz.ks = DOM_BUF;
    dom[i].Gfz.ksb = dom[i].Gfz.ks - DOM_BUF;
    dom[i].Gfz.kn = dom[i].zn + 1;
    dom[i].Gfz.knb = dom[i].Gfz.kn + 2 * DOM_BUF;
    dom[i].Gfz.ke = dom[i].Gfz.ks + dom[i].Gfz.kn;
    dom[i].Gfz.keb = dom[i].Gfz.ke + DOM_BUF;

    dom[i].Gfz.s1 = dom[i].Gfz.in;
    dom[i].Gfz.s2 = dom[i].Gfz.s1 * dom[i].Gfz.jn;
    dom[i].Gfz.s3 = dom[i].Gfz.s2 * dom[i].Gfz.kn;
    dom[i].Gfz.s1b = dom[i].Gfz.inb;
    dom[i].Gfz.s2b = dom[i].Gfz.s1b * dom[i].Gfz.jnb;
    dom[i].Gfz.s3b = dom[i].Gfz.s2b * dom[i].Gfz.knb;

    dom[i].Gfz._is = DOM_BUF;
    dom[i].Gfz._isb = dom[i].Gfz._is - DOM_BUF;
    dom[i].Gfz._in = dom[i].xn;
    dom[i].Gfz._inb = dom[i].Gfz._in + 2 * DOM_BUF;
    dom[i].Gfz._ie = dom[i].Gfz._is + dom[i].Gfz._in;
    dom[i].Gfz._ieb = dom[i].Gfz._ie + DOM_BUF;

    dom[i].Gfz._js = DOM_BUF;
    dom[i].Gfz._jsb = dom[i].Gfz._js - DOM_BUF;
    dom[i].Gfz._jn = dom[i].yn;
    dom[i].Gfz._jnb = dom[i].Gfz._jn + 2 * DOM_BUF;
    dom[i].Gfz._je = dom[i].Gfz._js + dom[i].Gfz._jn;
    dom[i].Gfz._jeb = dom[i].Gfz._je + DOM_BUF;

    dom[i].Gfz._ks = DOM_BUF;
    dom[i].Gfz._ksb = dom[i].Gfz._ks - DOM_BUF;
    dom[i].Gfz._kn = dom[i].zn + 1;
    dom[i].Gfz._knb = dom[i].Gfz._kn + 2 * DOM_BUF;
    dom[i].Gfz._ke = dom[i].Gfz._ks + dom[i].Gfz._kn;
    dom[i].Gfz._keb = dom[i].Gfz._ke + DOM_BUF;

    dom[i].Gfz._s1 = dom[i].Gfz._in;
    dom[i].Gfz._s2 = dom[i].Gfz._s1 * dom[i].Gfz._jn;
    dom[i].Gfz._s3 = dom[i].Gfz._s2 * dom[i].Gfz._kn;
    dom[i].Gfz._s1b = dom[i].Gfz._inb;
    dom[i].Gfz._s2b = dom[i].Gfz._s1b * dom[i].Gfz._jnb;
    dom[i].Gfz._s3b = dom[i].Gfz._s2b * dom[i].Gfz._knb;
  }

  // set up grid index structs
  // allocate and initialize pressure and velocity vectors
  p0 = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  p = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  divU = (real*) malloc(Dom.Gcc.s3b * sizeof(real));
  cpumem += Dom.Gcc.s3b * sizeof(real);
  u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u0 = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v0 = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w0 = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  diff0_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  diff0_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  diff0_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  conv0_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  conv0_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  conv0_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  diff_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  diff_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  diff_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  conv_u = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  conv_v = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  conv_w = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u_star = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  v_star = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  w_star = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);
  u_WE = (real*) malloc(Dom.Gfx.jnb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.jnb*Dom.Gfx.knb * sizeof(real);
  u_SN = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.knb * sizeof(real);
  u_BT = (real*) malloc(Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real));
  cpumem += Dom.Gfx.inb*Dom.Gfx.jnb * sizeof(real);
  v_WE = (real*) malloc(Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.jnb*Dom.Gfy.knb * sizeof(real);
  v_SN = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.knb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.knb * sizeof(real);
  v_BT = (real*) malloc(Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real));
  cpumem += Dom.Gfy.inb*Dom.Gfy.jnb * sizeof(real);
  w_WE = (real*) malloc(Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.jnb*Dom.Gfz.knb * sizeof(real);
  w_SN = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.knb * sizeof(real);
  w_BT = (real*) malloc(Dom.Gfz.inb*Dom.Gfz.jnb * sizeof(real));
  cpumem += Dom.Gfz.inb*Dom.Gfz.jnb * sizeof(real);
  f_x = (real*) malloc(Dom.Gfx.s3b * sizeof(real));
  cpumem += Dom.Gfx.s3b * sizeof(real);
  f_y = (real*) malloc(Dom.Gfy.s3b * sizeof(real));
  cpumem += Dom.Gfy.s3b * sizeof(real);
  f_z = (real*) malloc(Dom.Gfz.s3b * sizeof(real));
  cpumem += Dom.Gfz.s3b * sizeof(real);

  // set up the random number generator
  srand(time(NULL));

  for(i = 0; i < Dom.Gcc.s3b; i++) {
    p0[i] = 0.;
    p[i] = 0.;
    divU[i] = 0.;
  }
  for(i = 0; i < Dom.Gfx.s3b; i++) {
    u[i] = 0.;
    diff0_u[i] = 0.;
    diff_u[i] = 0.;
    u0[i] = 0.;
    conv0_u[i] = 0.;
    conv_u[i] = 0.;
    f_x[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.s3b; i++) {
    v[i] = 0.;
    diff0_v[i] = 0.;
    diff_v[i] = 0.;
    v0[i] = 0.;
    conv0_v[i] = 0.;
    conv_v[i] = 0.;
    f_y[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.s3b; i++) {
    w[i] = 0.;
    diff0_w[i] = 0.;
    diff_w[i] = 0.;
    w0[i] = 0.;
    conv0_w[i] = 0.;
    conv_w[i] = 0.;
    f_z[i] = 0.;
  }

  // integral scale
  turbl = (Dom.xl + Dom.yl + Dom.zl) / 3.;
  real urms = 3*turbA*turbl;

  // randomly initialize velocity components
  for(i = 0; i < Dom.Gfx.s3b; i++) {
    tmp = (rand() / (real)RAND_MAX - 0.5) * urms;
    u[i] = tmp;
    u0[i] = tmp;
  }
  for(i = 0; i < Dom.Gfy.s3b; i++) {
    tmp = (rand() / (real)RAND_MAX - 0.5) * urms;
    v[i] = tmp;
    v0[i] = tmp;
  }
  for(i = 0; i < Dom.Gfz.s3b; i++) {
    tmp = (rand() / (real)RAND_MAX - 0.5) * urms;
    w[i] = tmp;
    w0[i] = tmp;
  }

  // calculate the divergence of U
  real vol = (Dom.xn+2*DOM_BUF)*(Dom.yn+2*DOM_BUF)*(Dom.zn+2*DOM_BUF);
  for(k = Dom.Gcc.ksb; k < Dom.Gcc.keb; k++) {
    for(j = Dom.Gcc.jsb; j < Dom.Gcc.jeb; j++) {
      for(i = Dom.Gcc.isb; i < Dom.Gcc.ieb; i++) {
        W = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        E = (i+1) + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        S = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        N = i + (j+1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        B = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        T = i + j*Dom.Gfz.s1b + (k+1)*Dom.Gfz.s2b;
        C = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        p[C] = (u[E]-u[W])/Dom.dx + (v[N]-v[S])/Dom.dy + (w[T]-w[B])/Dom.dz;
        p[C] = p[C] / vol;
      }
    }
  }

  real umean = 0.;
  real vmean = 0.;
  real wmean = 0.;
  // subtract off the divergence of U to make U solenoidal
  for(k = Dom.Gfx.ks; k < Dom.Gfx.ke; k++) {
    for(j = Dom.Gfx.js; j < Dom.Gfx.je; j++) {
      for(i = Dom.Gfx.is; i < Dom.Gfx.ie; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        W = (i-1) + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        E = i + j*Dom.Gcc.s1b + k*Dom.Gcc.s2b;
        u[C] = u[C] - 0.5*(p[W] + p[E]);
        u0[C] = u0[C] - 0.5*(p[W] + p[E]);
        umean += u[C];
      }
    }
  }
  for(k = Dom.Gfy.ks; k < Dom.Gfy.ke; k++) {
    for(j = Dom.Gfy.js; j < Dom.Gfy.je; j++) {
      for(i = Dom.Gfy.is; i < Dom.Gfy.ie; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        S = i + (j-1)*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        N = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v[C] = v[C] - 0.5*(p[S] + p[N]);
        v0[C] = v0[C] - 0.5*(p[S] + p[N]);
        vmean += v[C];
      }
    }
  }
  for(k = Dom.Gfz.ks; k < Dom.Gfz.ke; k++) {
    for(j = Dom.Gfz.js; j < Dom.Gfz.je; j++) {
      for(i = Dom.Gfz.is; i < Dom.Gfz.ie; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        B = i + j*Dom.Gfz.s1b + (k-1)*Dom.Gfz.s2b;
        T = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w[C] = w[C] - 0.5*(p[B] + p[T]);
        w0[C] = w0[C] - 0.5*(p[B] + p[T]);
        wmean += w[C];
      }
    }
  }

  umean /= vol;
  vmean /= vol;
  wmean /= vol;

  // re-scale to give zero mean velocity in each direction
  for(k = Dom.Gfx.ksb; k < Dom.Gfx.keb; k++) {
    for(j = Dom.Gfx.jsb; j < Dom.Gfx.jeb; j++) {
      for(i = Dom.Gfx.isb; i < Dom.Gfx.ieb; i++) {
        C = i + j*Dom.Gfx.s1b + k*Dom.Gfx.s2b;
        u[C] = u[C] - umean;
        u0[C] = u0[C] - umean;
      }
    }
  }
  for(k = Dom.Gfy.ksb; k < Dom.Gfy.keb; k++) {
    for(j = Dom.Gfy.jsb; j < Dom.Gfy.jeb; j++) {
      for(i = Dom.Gfy.isb; i < Dom.Gfy.ieb; i++) {
        C = i + j*Dom.Gfy.s1b + k*Dom.Gfy.s2b;
        v[C] = v[C] - vmean;
        v0[C] = v0[C] - vmean;
      }
    }
  }
  for(k = Dom.Gfz.ksb; k < Dom.Gfz.keb; k++) {
    for(j = Dom.Gfz.jsb; j < Dom.Gfz.jeb; j++) {
      for(i = Dom.Gfz.isb; i < Dom.Gfz.ieb; i++) {
        C = i + j*Dom.Gfz.s1b + k*Dom.Gfz.s2b;
        w[C] = w[C] - wmean;
        w0[C] = w0[C] - wmean;
      }
    }
  }

  for(i = 0; i < Dom.Gfx.jnb*Dom.Gfx.knb; i++) {
    u_WE[i] = 0.;
  }
  for(i = 0; i < Dom.Gfx.inb*Dom.Gfx.knb; i++) {
    u_SN[i] = 0.;
  }
  for(i = 0; i < Dom.Gfx.inb*Dom.Gfx.jnb; i++) {
    u_BT[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.jnb*Dom.Gfy.knb; i++) {
    v_WE[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.inb*Dom.Gfy.knb; i++) {
    v_SN[i] = 0.;
  }
  for(i = 0; i < Dom.Gfy.inb*Dom.Gfy.jnb; i++) {
    v_BT[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.jnb*Dom.Gfz.knb; i++) {
    w_WE[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.inb*Dom.Gfz.knb; i++) {
    w_SN[i] = 0.;
  }
  for(i = 0; i < Dom.Gfz.inb*Dom.Gfz.jnb; i++) {
    w_BT[i] = 0.;
  }

  // initialize some variables
  dt = 2 * nu / (Dom.dx * Dom.dx);
  dt += 2 * nu / (Dom.dy * Dom.dy);
  dt += 2 * nu / (Dom.dz * Dom.dz);
  dt = CFL / dt;
  dt0 = -1.;
  stepnum = 0;
  rec_flow_field_stepnum_out = 0;
  rec_paraview_stepnum_out = 0;
  rec_point_particle_stepnum_out = 0;

  return EXIT_SUCCESS;
}

