// Nutthapat Sumanont, personal site (geometry version).
// Story: 276 shapes, one per patient in the Beat AML 2.0 cohort used in the NPM1 project
// (100 NPM1-mutated, 176 wild-type; 37 of the 100 in the primitive-like subgroup). 111 more
// drop in for the TCGA-LAML replication (31 NPM1-mutated, 80 wild-type). Only group sizes
// are used, never patient data. The sesamin panel uses Beck's mean MTT viability values.
(function () {
  "use strict";

  var PI = Math.PI, TAU = PI * 2, SQ3 = Math.sqrt(3);
  var reduce = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;

  var seed = 2026;
  function rand() {
    seed |= 0; seed = (seed + 0x6D2B79F5) | 0;
    var t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }
  function clamp01(v) { return v < 0 ? 0 : v > 1 ? 1 : v; }
  function mix(a, b, t) { return a + (b - a) * t; }
  function ease(t) { return t < 0.5 ? 4 * t * t * t : 1 - Math.pow(-2 * t + 2, 3) / 2; }
  function bounce(t) {
    var n = 7.5625, d = 2.75;
    if (t < 1 / d) return n * t * t;
    if (t < 2 / d) return n * (t -= 1.5 / d) * t + 0.75;
    if (t < 2.5 / d) return n * (t -= 2.25 / d) * t + 0.9375;
    return n * (t -= 2.625 / d) * t + 0.984375;
  }
  function mixc(a, b, t) { return [mix(a[0], b[0], t), mix(a[1], b[1], t), mix(a[2], b[2], t)]; }
  function css(c, a) { return "rgba(" + (c[0] | 0) + "," + (c[1] | 0) + "," + (c[2] | 0) + "," + a.toFixed(3) + ")"; }
  function tokenRGB(name, fallback) {
    var v = getComputedStyle(document.documentElement).getPropertyValue(name).trim() || fallback;
    var m = v.match(/^#([0-9a-f]{6})$/i);
    if (!m) return [128, 128, 128];
    var n = parseInt(m[1], 16); return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
  }
  var FONT = "-apple-system, BlinkMacSystemFont, 'Helvetica Neue', 'Segoe UI', Roboto, Arial, sans-serif";

  // Round clusters of pointy-top hexagons (unit spacing, hexagon radius 1), closest to the centre first
  function hexCluster(n) {
    var pts = [];
    for (var a = -15; a <= 15; a++) for (var b = -15; b <= 15; b++) {
      var x = SQ3 * (a + b / 2), y = 1.5 * b;
      pts.push([x, y, Math.round(Math.hypot(x, y) * 1000), Math.atan2(y, x)]);
    }
    pts.sort(function (p1, p2) { return p1[2] - p2[2] || p1[3] - p2[3]; });
    var out = pts.slice(0, n).map(function (p) { return [p[0], p[1]]; });
    var x0 = 1e9, x1 = -1e9, y0 = 1e9, y1 = -1e9;
    out.forEach(function (p) { x0 = Math.min(x0, p[0] - 1); x1 = Math.max(x1, p[0] + 1); y0 = Math.min(y0, p[1] - 1); y1 = Math.max(y1, p[1] + 1); });
    return { pts: out, w: x1 - x0, h: y1 - y0, cx: (x0 + x1) / 2, cy: (y0 + y1) / 2 };
  }

  /* =====================================================================
     1. The story canvas
     ===================================================================== */
  (function story() {
    var canvas = document.getElementById("cells");
    var storyEl = document.getElementById("story");
    if (!canvas || !canvas.getContext) return;
    var ctx = canvas.getContext("2d");
    var steps = Array.prototype.slice.call(document.querySelectorAll(".story-steps .step"));
    var LAST = steps.length - 1;
    var heroText = document.querySelector(".step-hero .step-text");
    var card = document.getElementById("project");

    var N_MAIN = 276, N_MUT = 100, N_PRIM = 37, N_T = 111, N_TMUT = 31;
    var SAMPLES = 24;
    function profile(m, phi, scale) {
      var out = new Float32Array(SAMPLES), seg = TAU / m;
      for (var k = 0; k < SAMPLES; k++) {
        var th = -PI / 2 + TAU * k / SAMPLES, a = (((th - phi) % seg) + seg) % seg;
        out[k] = scale * Math.cos(PI / m) / Math.cos(a - PI / m);
      }
      return out;
    }
    var HEX = profile(6, -PI / 2, 1), SQUARE = profile(4, PI / 4, 1);
    var KINDS = [profile(3, -PI / 2, 1.3), profile(4, PI / 4, 1.12), profile(5, -PI / 2, 1.04), HEX];
    var COS = new Float32Array(SAMPLES), SIN = new Float32Array(SAMPLES);
    for (var q = 0; q < SAMPLES; q++) { var th = -PI / 2 + TAU * q / SAMPLES; COS[q] = Math.cos(th); SIN[q] = Math.sin(th); }
    var mixed = new Float32Array(SAMPLES);
    function mixProfile(a, b, m) { for (var k = 0; k < SAMPLES; k++) mixed[k] = mix(a[k], b[k], m); return mixed; }

    var C_ALL = hexCluster(N_MAIN), C_MUT = hexCluster(N_MUT), C_WT = hexCluster(N_MAIN - N_MUT);
    var C_PRIM = hexCluster(N_PRIM), C_MAT = hexCluster(N_MUT - N_PRIM);
    var C_TMUT = hexCluster(N_TMUT), C_TWT = hexCluster(N_T - N_TMUT);

    var cells = [], i;
    for (i = 0; i < N_MAIN; i++) cells.push(makeCell(0, i < N_MUT, i < N_MUT ? i : i - N_MUT));
    for (i = 0; i < N_T; i++) cells.push(makeCell(1, i < N_TMUT, i < N_TMUT ? i : i - N_TMUT));
    function makeCell(cohort, mut, rank) {
      var sub = mut && cohort === 0 ? (rank < N_PRIM ? 0 : 1) : -1;
      return {
        cohort: cohort, mut: mut, rank: rank, sub: sub, subRank: sub === 1 ? rank - N_PRIM : rank,
        prof: KINDS[Math.floor(rand() * KINDS.length)], outline: rand() < 0.45, tone: 0.18 + rand() * 0.62,
        fx: 0, fy: 0, vx: (rand() - 0.5) * 0.8, vy: (rand() - 0.5) * 0.8,
        rot: rand() * TAU, vrot: (rand() - 0.5) * 0.02, size: 0.7 + rand() * 0.9,
        x: 0, y: 0, a: 0, wd: 0, wd2: 0, delay: rand() * 0.45, t: []
      };
    }
    var order = [];
    for (i = 0; i < N_MAIN; i++) order.push(i);
    for (i = N_MAIN - 1; i > 0; i--) { var j = Math.floor(rand() * (i + 1)); var tmp = order[i]; order[i] = order[j]; order[j] = tmp; }
    var slot = []; for (i = 0; i < N_MAIN; i++) slot[order[i]] = i;

    var col = {};
    function readColours() {
      col.cell = tokenRGB("--cell", "#1d1d1f"); col.grey = tokenRGB("--cell-grey", "#c7c7cc");
      col.orange = tokenRGB("--cell-orange", "#ff6a00"); col.orange2 = tokenRGB("--cell-orange-2", "#ffc29e");
      col.glow = tokenRGB("--cell-glow", "#ffb36b"); col.label = tokenRGB("--story-ink-2", "#6e6e73");
      col.ink = tokenRGB("--story-ink", "#1d1d1f");
    }

    var W = 0, H = 0, dpr = 1, area = null, R0 = 8, RS = [], labels = {};

    function heroRect() {
      var r = { left: 1e9, right: -1e9, top: 1e9, bottom: -1e9 };
      Array.prototype.forEach.call(heroText.children, function (el) {
        var b = el.getBoundingClientRect();
        if (el.tagName === "H1") { var range = document.createRange(); range.selectNodeContents(el); b = range.getBoundingClientRect(); }
        r.left = Math.min(r.left, b.left); r.right = Math.max(r.right, b.right); r.top = Math.min(r.top, b.top); r.bottom = Math.max(r.bottom, b.bottom);
      });
      return r;
    }

    function layout() {
      var rect = canvas.getBoundingClientRect();
      W = rect.width; H = rect.height;
      dpr = Math.min(window.devicePixelRatio || 1, 2);
      canvas.width = Math.round(W * dpr); canvas.height = Math.round(H * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

      var wide = W >= 900;
      area = wide ? { x: W * 0.45, y: 72, w: W * 0.5, h: H - 124 } : { x: 16, y: 64, w: W - 32, h: H * 0.58 - 64 };
      var cx = area.x + area.w / 2, cy = area.y + area.h / 2;
      R0 = Math.max(6, Math.min(12, Math.min(W, H) * 0.012));
      function place(cl, k, ox, oy, sc) { var p = cl.pts[k]; return [ox + (p[0] - cl.cx) * sc, oy + (p[1] - cl.cy) * sc]; }
      var gap = 3, lab = wide ? 58 : 68;

      // States 1 to 4: one cluster, then NPM1-mutated next to wild-type
      var pairW = C_MUT.w + gap + C_WT.w, pairH = Math.max(C_MUT.h, C_WT.h);
      var s1 = Math.min(area.w / Math.max(C_ALL.w, pairW), (area.h - 56) / Math.max(C_ALL.h, pairH), 13);
      var yPair = cy - 14;
      var xM = cx - pairW * s1 / 2 + C_MUT.w * s1 / 2, xW = xM + (C_MUT.w / 2 + gap + C_WT.w / 2) * s1;

      // State 5: Beat AML on top, TCGA below
      var rowGap = 3.4, botW = C_TMUT.w + gap + C_TWT.w, botH = Math.max(C_TMUT.h, C_TWT.h);
      var s5 = Math.min(area.w / Math.max(pairW, botW), (area.h - 64) / (pairH + rowGap + botH), 10);
      var top5 = cy - ((pairH + rowGap + botH) * s5 + 40) / 2 + 18;
      var yTop = top5 + pairH * s5 / 2, yBot = top5 + (pairH + rowGap) * s5 + 22 + botH * s5 / 2;
      var xM5 = cx - pairW * s5 / 2 + C_MUT.w * s5 / 2, xW5 = xM5 + (C_MUT.w / 2 + gap + C_WT.w / 2) * s5;
      var xTM = cx - botW * s5 / 2 + C_TMUT.w * s5 / 2, xTW = xTM + (C_TMUT.w / 2 + gap + C_TWT.w / 2) * s5;

      // States 6 to 8: the 100 split into 37 primitive-like and 63 more mature, next to wild-type.
      // On a wide screen the last step slides the clusters right, clear of the project card.
      function trio(ax, aw) {
        var o = {}, acx = ax + aw / 2;
        if (wide) {
          var trioW = C_PRIM.w + C_MAT.w + C_WT.w + 2 * gap;
          o.s = Math.min(aw / trioW, (area.h - lab - 20) / C_WT.h, 12);
          var x0 = acx - trioW * o.s / 2;
          o.xP = x0 + C_PRIM.w * o.s / 2; o.xQ = x0 + (C_PRIM.w + gap) * o.s + C_MAT.w * o.s / 2;
          o.xR = x0 + (C_PRIM.w + gap + C_MAT.w + gap) * o.s + C_WT.w * o.s / 2;
          o.yP = o.yQ = o.yR = cy - lab / 2;
        } else {
          var L = Math.max(C_PRIM.w, C_MAT.w), twoW = L + gap + C_WT.w;
          o.s = Math.min(aw / twoW, (area.h - 2 * lab - 6) / (C_PRIM.h + C_MAT.h), (area.h - lab) / C_WT.h, 12);
          o.xP = o.xQ = acx - twoW * o.s / 2 + L * o.s / 2; o.xR = o.xP + (L / 2 + gap + C_WT.w / 2) * o.s;
          var colTop = cy - ((C_PRIM.h + C_MAT.h) * o.s + 2 * lab) / 2;
          o.yP = colTop + C_PRIM.h * o.s / 2; o.yQ = colTop + C_PRIM.h * o.s + lab + C_MAT.h * o.s / 2; o.yR = cy - lab / 2;
        }
        return o;
      }
      var T6 = trio(area.x, area.w), T8 = T6;
      if (wide && card) {
        var cardRight = card.getBoundingClientRect().right - rect.left;
        var x8 = Math.max(W * 0.58, cardRight + 40);
        T8 = trio(x8, W - x8 - 64);
      }
      RS = [s1, s1, s1, s1, s1, s5, T6.s, T6.s, T8.s].map(function (s) { return s * 0.9; });

      cells.forEach(function (c, idx) {
        var p1, p2, p5, p6, p8;
        if (c.cohort === 0) {
          p1 = place(C_ALL, slot[idx], cx, cy, s1);
          if (c.mut) {
            p2 = place(C_MUT, c.rank, xM, yPair, s1); p5 = place(C_MUT, c.rank, xM5, yTop, s5);
            p6 = c.sub === 0 ? place(C_PRIM, c.subRank, T6.xP, T6.yP, T6.s) : place(C_MAT, c.subRank, T6.xQ, T6.yQ, T6.s);
            p8 = c.sub === 0 ? place(C_PRIM, c.subRank, T8.xP, T8.yP, T8.s) : place(C_MAT, c.subRank, T8.xQ, T8.yQ, T8.s);
            c.wd = Math.hypot(C_MUT.pts[c.rank][0], C_MUT.pts[c.rank][1]);
            var sp = (c.sub === 0 ? C_PRIM : C_MAT).pts[c.subRank]; c.wd2 = Math.hypot(sp[0], sp[1]);
          } else {
            p2 = place(C_WT, c.rank, xW, yPair, s1); p5 = place(C_WT, c.rank, xW5, yTop, s5);
            p6 = place(C_WT, c.rank, T6.xR, T6.yR, T6.s); p8 = place(C_WT, c.rank, T8.xR, T8.yR, T8.s);
          }
        } else {
          p5 = c.mut ? place(C_TMUT, c.rank, xTM, yBot, s5) : place(C_TWT, c.rank, xTW, yBot, s5);
          if (c.mut) c.wd = Math.hypot(C_TMUT.pts[c.rank][0], C_TMUT.pts[c.rank][1]);
          p1 = p2 = p6 = p8 = p5;
        }
        c.t = [null, p1, p2, p2, p2, p5, p6, p6, p8];
        if (!c.fx && !c.fy) { c.fx = 20 + rand() * (W - 40); c.fy = 20 + rand() * (H - 40); c.x = c.fx; c.y = c.fy; c.a = c.rot; }
      });

      labels.pair = [
        { x: xM, y: yPair + C_MUT.h * s1 / 2 + 26, text: "NPM1-mutated", n: "100", orange: true },
        { x: xW, y: yPair + C_WT.h * s1 / 2 + 26, text: "Wild-type", n: "176" }
      ];
      labels.rows = [
        { x: cx, y: yTop - pairH * s5 / 2 - 12, text: "Beat AML 2.0, 276 patients" },
        { x: cx, y: yBot - botH * s5 / 2 - 12, text: "TCGA-LAML, 111 patients" }
      ];
      function trioLabels(T) {
        return [
          { x: T.xP, y: T.yP + C_PRIM.h * T.s / 2 + 22, text: "More primitive", n: "37", orange: true, sub: wide ? "more sensitive to venetoclax" : "more sensitive" },
          { x: T.xQ, y: T.yQ + C_MAT.h * T.s / 2 + 22, text: "More mature", n: "63", orange: true, sub: "less sensitive" },
          { x: T.xR, y: T.yR + C_WT.h * T.s / 2 + 22, text: "Wild-type", n: "176" }
        ];
      }
      labels.trio6 = trioLabels(T6); labels.trio8 = trioLabels(T8);
      egg.targets = null;
    }

    function progress() {
      var tops = steps.map(function (el) { return el.getBoundingClientRect().top; });
      if (tops[0] >= 0) return 0;
      for (var k = 0; k < LAST; k++) {
        if (tops[k + 1] > 0) {
          var span = tops[k + 1] - tops[k];
          return k + Math.min(1, Math.max(0, -tops[k] / (span || window.innerHeight)));
        }
      }
      return LAST;
    }

    // Play: stir with the cursor, tap to scatter, double-tap to spell BECK
    var px = -9999, py = -9999, pvx = 0, pvy = 0, lastMove = 0, burst = 0, bx = 0, by = 0, lastTap = { t: 0, x: 0, y: 0 };
    var egg = { start: -1, targets: null, size: 3 };
    function toLocal(e) { var r = canvas.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; }
    window.addEventListener("pointermove", function (e) {
      var p = toLocal(e), now = performance.now(), dt = Math.max(8, now - lastMove);
      if (px > -9000) { pvx = (p[0] - px) / dt * 16; pvy = (p[1] - py) / dt * 16; }
      px = p[0]; py = p[1]; lastMove = now;
    }, { passive: true });
    storyEl.addEventListener("pointerdown", function (e) {
      if (e.target.closest && e.target.closest("a, button")) return;
      var p = toLocal(e), now = performance.now();
      bx = p[0]; by = p[1]; burst = 1;
      if (now - lastTap.t < 330 && Math.hypot(p[0] - lastTap.x, p[1] - lastTap.y) < 50 && progress() < 0.5) startEgg();
      lastTap = { t: now, x: p[0], y: p[1] };
    }, { passive: true });
    document.addEventListener("pointerup", function (e) { if (e.pointerType !== "mouse") { px = py = -9999; pvx = pvy = 0; } }, { passive: true });

    function startEgg() { if (reduce) return; egg.targets = wordTargets("BECK"); egg.start = performance.now(); }
    function wordTargets(word) {
      var box = { x: area.x, y: area.y, w: area.w, h: area.h };
      if (heroText) {
        var cr = canvas.getBoundingClientRect(), tr = heroRect();
        if (W >= 900) { box.x = Math.max(area.x, tr.right - cr.left + 32); box.w = W - box.x - 40; }
        else { box.y = 64; box.h = Math.max(120, tr.top - cr.top - 80); }
      }
      var off = document.createElement("canvas");
      var w = Math.max(40, Math.round(box.w)), h = Math.max(40, Math.round(Math.min(box.h, box.w * 0.6)));
      off.width = w; off.height = h;
      var g = off.getContext("2d");
      var fs = Math.min(w / (word.length * 0.72), h * 0.9);
      g.font = "800 " + fs + "px " + FONT; g.textAlign = "center"; g.textBaseline = "middle"; g.fillStyle = "#000";
      g.fillText(word, w / 2, h / 2);
      var data = g.getImageData(0, 0, w, h).data, pts = [], step = 4;
      for (var tries = 0; tries < 12; tries++) {
        pts = [];
        for (var y = 0; y < h; y += step) for (var x = 0; x < w; x += step) if (data[(y * w + x) * 4 + 3] > 128) pts.push([x, y]);
        if (pts.length < N_MAIN * 1.15) break;
        step += 1;
      }
      var ox = box.x + box.w / 2 - w / 2, oy = box.y + box.h / 2 - h / 2, out = [];
      for (var k = 0; k < N_MAIN; k++) { var p = pts[Math.floor(k * pts.length / N_MAIN)] || [w / 2, h / 2]; out.push([ox + p[0], oy + p[1]]); }
      egg.size = step * 0.5;
      return out;
    }
    function eggAmount(now) {
      if (egg.start < 0) return 0;
      var t = (now - egg.start) / 1000;
      if (t < 0.8) return ease(t / 0.8);
      if (t < 3.2) return 1;
      if (t < 4.1) return 1 - ease((t - 3.2) / 0.9);
      egg.start = -1; return 0;
    }

    // Pager: one hexagon per step (wide screens)
    var pager = document.createElement("nav");
    pager.className = "pager"; pager.setAttribute("aria-label", "Story steps");
    var pagerButtons = steps.map(function (el, k) {
      var b = document.createElement("button");
      b.type = "button"; b.setAttribute("aria-label", "Go to part " + (k + 1) + " of " + steps.length);
      b.innerHTML = '<svg viewBox="0 0 12 14" aria-hidden="true"><polygon points="6,0 12,3.5 12,10.5 6,14 0,10.5 0,3.5"/></svg>';
      b.addEventListener("click", function () { window.scrollTo({ top: el.getBoundingClientRect().top + window.scrollY + 2, behavior: reduce ? "auto" : "smooth" }); });
      pager.appendChild(b); return b;
    });
    document.body.appendChild(pager);
    var pagerActive = -1;

    function shapePath(x, y, rad, ang, sx, prof, morphTo, m) {
      var ca = Math.cos(ang), sa = Math.sin(ang);
      ctx.beginPath();
      for (var k = 0; k < SAMPLES; k++) {
        var r = rad * mix(prof[k], morphTo[k], m), lx = COS[k] * r * sx, ly = SIN[k] * r;
        var X = x + lx * ca - ly * sa, Y = y + lx * sa + ly * ca;
        if (k === 0) ctx.moveTo(X, Y); else ctx.lineTo(X, Y);
      }
      ctx.closePath();
    }

    var running = false, last = 0, time = 0;
    function frame(now) {
      if (!running) return;
      var dt = Math.min(48, now - (last || now)); last = now; time += dt / 1000;
      draw(dt, now);
      requestAnimationFrame(frame);
    }

    function draw(dt, now) {
      var s = progress();
      ctx.clearRect(0, 0, W, H);
      var k = Math.min(LAST - 1, Math.floor(s)), t = ease(s - k);
      if (s >= LAST) { k = LAST - 1; t = 1; }
      var f60 = dt / 16.7;
      var spring = reduce ? 1 : 1 - Math.pow(0.84, f60);
      var morph = ease(clamp01(s));
      var cellR = s < 1 ? RS[1] : mix(RS[Math.max(1, k)], RS[Math.max(1, k + 1)], t);
      var greyT = clamp01((s - 1.3) / 0.5);
      var wave = clamp01((s - 2.35) / 0.5) * (1 - 0.65 * clamp01((s - 3.5) / 0.6)) * (1 - clamp01((s - 4.5) / 0.4));
      var tcgaWave = 0.35 * clamp01((s - 4.9) / 0.4);
      var tcgaOut = clamp01((s - 5.25) / 0.5);
      var split = clamp01((s - 5.35) / 0.6);
      var venWave = clamp01((s - 6.35) / 0.5);
      var wtFade = mix(1, 0.5, split);
      var dim = W < 900 ? 1 - 0.75 * clamp01((s - (LAST - 0.7)) / 0.5) : 1;
      var eg = s < 0.6 ? eggAmount(now) : (egg.start = -1, 0);
      burst *= reduce ? 0 : Math.pow(0.9, f60);
      pvx *= Math.pow(0.85, f60); pvy *= Math.pow(0.85, f60);

      var hr = null;
      if (s < 1 && heroText) {
        var cr = canvas.getBoundingClientRect(), tr = heroRect();
        hr = { l: tr.left - cr.left - 12, r: tr.right - cr.left + 12, t: tr.top - cr.top - 12, b: tr.bottom - cr.top + 12 };
      }

      for (var n = 0; n < cells.length; n++) {
        var c = cells[n];

        if (!reduce) {   // hero physics: drift, stir, scatter, bounce off the edges
          c.vx += (rand() - 0.5) * 0.04 * f60; c.vy += (rand() - 0.5) * 0.04 * f60;
          var ddx = c.fx - px, ddy = c.fy - py, dd = Math.hypot(ddx, ddy);
          if (s < 1 && dd < 120) {
            var w = 1 - dd / 120;
            c.vx += pvx * 0.06 * w + (ddx / (dd || 1)) * 0.25 * w; c.vy += pvy * 0.06 * w + (ddy / (dd || 1)) * 0.25 * w;
            c.vrot += (pvx - pvy) * 0.0006 * w;
          }
          if (burst > 0.5 && s < 1) {
            var ex = c.fx - bx, ey = c.fy - by, ed = Math.hypot(ex, ey) || 1;
            if (ed < 220) { var bf = (1 - ed / 220) * 1.6 * burst; c.vx += ex / ed * bf; c.vy += ey / ed * bf; c.vrot += (rand() - 0.5) * 0.06 * burst; }
          }
          c.vx *= Math.pow(0.975, f60); c.vy *= Math.pow(0.975, f60); c.vrot *= Math.pow(0.99, f60);
          var spd = Math.hypot(c.vx, c.vy); if (spd > 9) { c.vx *= 9 / spd; c.vy *= 9 / spd; }
          c.fx += c.vx * f60; c.fy += c.vy * f60; c.rot += (c.vrot + 0.002) * f60;
          var m0 = R0 * 1.6;
          if (c.fx < m0) { c.fx = m0; c.vx = Math.abs(c.vx); } if (c.fx > W - m0) { c.fx = W - m0; c.vx = -Math.abs(c.vx); }
          if (c.fy < m0) { c.fy = m0; c.vy = Math.abs(c.vy); } if (c.fy > H - m0) { c.fy = H - m0; c.vy = -Math.abs(c.vy); }
        }

        var tx, ty;
        if (s < 1 && c.cohort === 0) { tx = mix(c.fx, c.t[1][0], morph); ty = mix(c.fy, c.t[1][1], morph); }
        else { var A = c.t[Math.max(1, k)], B = c.t[Math.max(1, k + 1)]; tx = mix(A[0], B[0], t); ty = mix(A[1], B[1], t); }

        var drop = 1, alpha = 1;
        if (c.cohort === 1) {   // TCGA drops in, then falls away again for the subgroup steps
          drop = clamp01((s - 4.15 - c.delay * 0.9) / 0.45);
          if (drop <= 0 || tcgaOut >= 1) continue;
          ty = ty - (1 - bounce(drop)) * H * 0.55 + ease(tcgaOut) * H * 0.35;
          alpha = 1 - tcgaOut;
        }
        if (c.cohort === 0 && eg > 0 && egg.targets) { var E = egg.targets[slot[n]]; tx = mix(tx, E[0], eg); ty = mix(ty, E[1], eg); }

        if (s >= 1) {
          var pd = Math.hypot(tx - px, ty - py);
          if (pd < 70 && pd > 0.01) { var pf = (1 - pd / 70) * 18; tx += (tx - px) / pd * pf; ty += (ty - py) / pd * pf; }
          if (burst > 0.01) { var gx = tx - bx, gy = ty - by, gd = Math.hypot(gx, gy) || 1; if (gd < 240) { var gf = (1 - gd / 240) * 60 * burst; tx += gx / gd * gf; ty += gy / gd * gf; } }
        }
        if (c.cohort === 1 && (drop < 1 || tcgaOut > 0)) { c.x = tx; c.y = ty; } else { c.x += (tx - c.x) * spring; c.y += (ty - c.y) * spring; }

        var target = mix(mix(c.rot, 0, morph), 0, c.cohort === 0 ? eg : 0);
        if (c.cohort === 0) { var da = Math.atan2(Math.sin(target - c.a), Math.cos(target - c.a)); c.a += da * (reduce ? 1 : Math.min(1, 0.2 * f60)); } else c.a = 0;

        var rad, profM = 1, sx = 1, colour, fillA = 1, strokeA = 0;
        if (c.cohort === 0) {
          rad = mix(R0 * c.size, cellR, morph); profM = morph;
          fillA = c.outline ? morph : 1; strokeA = c.outline ? 1 - morph : 0;
          colour = col.cell;
          if (c.mut) {
            var f = clamp01((s - 1.25 - c.delay) / 0.35);   // each cell flips over to orange
            sx = Math.abs(Math.cos(PI * f));
            if (f >= 0.5) {
              colour = c.sub === 1 ? mixc(col.orange, col.orange2, split) : col.orange;
              var amp = Math.max(wave, c.sub === 0 ? venWave : 0);
              if (amp > 0.01) {
                var pulse = 0.5 + 0.5 * Math.sin(time * 3.4 - (s > 5.5 ? c.wd2 : c.wd) * 0.85);
                colour = mixc(colour, col.glow, pulse * amp * 0.7);
                rad *= 1 + 0.14 * pulse * amp;
              }
            }
          } else {
            colour = mixc(col.cell, col.grey, greyT);
            alpha *= wtFade;
          }
          if (s < 1) alpha *= mix(c.tone, 1, morph);
        } else {
          rad = RS[5] * (0.6 + 0.4 * drop);
          colour = c.mut ? col.orange : col.grey;
          if (c.mut && tcgaWave > 0) { var p2 = 0.5 + 0.5 * Math.sin(time * 3.4 - c.wd * 0.85); colour = mixc(col.orange, col.glow, p2 * tcgaWave); rad *= 1 + 0.1 * p2 * tcgaWave; }
        }
        var prof = c.cohort === 0 ? c.prof : HEX, to = HEX;
        if (eg > 0 && c.cohort === 0) {
          to = SQUARE; profM = eg; prof = mixProfile(c.prof, HEX, morph);
          rad = mix(rad, egg.size, eg); colour = mixc(colour, col.orange, eg); fillA = mix(fillA, 1, eg); strokeA *= 1 - eg; alpha = mix(alpha, 1, eg);
        }
        alpha *= dim;
        if (hr && c.x > hr.l && c.x < hr.r && c.y > hr.t && c.y < hr.b) alpha *= mix(0.18, 1, clamp01(s * 1.5 + eg));

        shapePath(c.x, c.y, rad, c.a, sx, prof, to, profM);
        if (fillA > 0.01) { ctx.fillStyle = css(colour, alpha * fillA); ctx.fill(); }
        if (strokeA > 0.01) { ctx.lineWidth = 1.6; ctx.strokeStyle = css(colour, alpha * strokeA); ctx.stroke(); }
      }

      // Labels
      var pairA = clamp01((s - 1.6) / 0.4) * (1 - clamp01((s - 4.1) / 0.4)) * dim;
      var rowA = clamp01((s - 4.5) / 0.4) * (1 - clamp01((s - 5.2) / 0.3)) * dim;
      var trioA = clamp01((s - 5.6) / 0.4) * dim, venA = clamp01((s - 6.45) / 0.4) * dim;
      ctx.textAlign = "center";
      function group(L, a) {
        ctx.font = "600 13px " + FONT; ctx.fillStyle = css(L.orange ? col.orange : col.label, a); ctx.fillText(L.text, L.x, L.y);
        ctx.font = "700 20px " + FONT; ctx.fillStyle = css(col.ink, a); ctx.fillText(L.n, L.x, L.y + 24);
        if (L.sub && venA > 0) { ctx.font = "500 12px " + FONT; ctx.fillStyle = css(L.text === "More primitive" ? col.orange : col.label, venA); ctx.fillText(L.sub, L.x, L.y + 42); }
      }
      if (pairA > 0) labels.pair.forEach(function (L) { group(L, pairA); });
      if (trioA > 0) {
        var u = ease(clamp01(s - (LAST - 1)));
        labels.trio6.forEach(function (L, x) {
          var M = labels.trio8[x];
          group({ x: mix(L.x, M.x, u), y: mix(L.y, M.y, u), text: L.text, n: L.n, orange: L.orange, sub: L.sub }, trioA);
        });
      }
      if (rowA > 0) { ctx.font = "600 13px " + FONT; ctx.fillStyle = css(col.label, rowA); labels.rows.forEach(function (L) { ctx.fillText(L.text, L.x, L.y); }); }

      // Pager
      var rect = storyEl.getBoundingClientRect();
      pager.classList.toggle("show", rect.top < window.innerHeight * 0.5 && rect.bottom > window.innerHeight * 0.5);
      var act = Math.min(LAST, Math.round(s));
      if (act !== pagerActive) {
        pagerButtons.forEach(function (b, x) { if (x === act) b.setAttribute("aria-current", "step"); else b.removeAttribute("aria-current"); });
        pagerActive = act;
      }
    }

    function start() { if (!running) { running = true; last = 0; requestAnimationFrame(frame); } }
    function stop() { running = false; pager.classList.remove("show"); }

    readColours();
    layout();
    draw(16, performance.now());
    if ("IntersectionObserver" in window) {
      new IntersectionObserver(function (es) { es.forEach(function (en) { if (en.isIntersecting) start(); else stop(); }); }).observe(storyEl);
    } else start();
    document.addEventListener("visibilitychange", function () { if (document.hidden) stop(); else start(); });
    var mq = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;
    if (mq && mq.addEventListener) mq.addEventListener("change", readColours);
    new MutationObserver(readColours).observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
    var rt;
    window.addEventListener("resize", function () { clearTimeout(rt); rt = setTimeout(function () { layout(); if (!running) draw(16, performance.now()); }, 120); });
  })();

  /* =====================================================================
     2. Sesamin: two hexagon "wells" driven by the real mean viability values
     ===================================================================== */
  (function lab() {
    var root = document.getElementById("lab");
    if (!root) return;
    var DATA = {
      "NB4":    { 24: [100, 86.6, 75.3, 68.9], 48: [100, 77.7, 68.2, 53.7] },
      "MOLT-4": { 24: [100, 85.9, 70.8, 56.9], 48: [100, 73.9, 65.1, 52.4] }
    };
    var N = 60, state = { dose: 3, time: 48 };
    var cl = hexCluster(N), sc = 180 / Math.max(cl.w, cl.h), r = sc * 0.86;
    var NS = "http://www.w3.org/2000/svg";
    var wells = Array.prototype.slice.call(root.querySelectorAll(".well")).map(function (fig, wi) {
      var svg = fig.querySelector("svg"), polys = [];
      var perm = []; for (var k = 0; k < N; k++) perm.push(k);
      for (k = N - 1; k > 0; k--) { var j = Math.floor(rand() * (k + 1)); var tmp = perm[k]; perm[k] = perm[j]; perm[j] = tmp; }
      cl.pts.forEach(function (p, k) {
        var x = (p[0] - cl.cx) * sc, y = (p[1] - cl.cy) * sc, pts = [];
        for (var v = 0; v < 6; v++) { var a = -PI / 2 + v * PI / 3; pts.push((x + r * Math.cos(a)).toFixed(2) + "," + (y + r * Math.sin(a)).toFixed(2)); }
        var poly = document.createElementNS(NS, "polygon");
        poly.setAttribute("points", pts.join(" "));
        svg.appendChild(poly); polys.push(poly);
      });
      return { line: fig.getAttribute("data-line"), fig: fig, polys: polys, deathOrder: perm };
    });

    function render() {
      wells.forEach(function (w) {
        var v = DATA[w.line][state.time][state.dose];
        var dead = N - Math.round(v / 100 * N);
        w.deathOrder.forEach(function (idx, rank) {
          var p = w.polys[idx], isDead = rank < dead;
          p.style.transitionDelay = reduce ? "0ms" : (rank % 20) * 18 + "ms";
          p.classList.toggle("dead", isDead);
        });
        w.fig.querySelector(".pct").textContent = v.toFixed(1).replace(/\.0$/, "") + "%";
      });
      root.querySelectorAll("[data-dose]").forEach(function (b) { b.setAttribute("aria-pressed", String(+b.getAttribute("data-dose") === state.dose)); });
      root.querySelectorAll("[data-time]").forEach(function (b) { b.setAttribute("aria-pressed", String(+b.getAttribute("data-time") === state.time)); });
    }
    root.addEventListener("click", function (e) {
      var b = e.target.closest("button");
      if (!b) return;
      if (b.hasAttribute("data-dose")) state.dose = +b.getAttribute("data-dose");
      else if (b.hasAttribute("data-time")) state.time = +b.getAttribute("data-time");
      else return;
      render();
    });
    render();
  })();

  /* =====================================================================
     3. Contact: a honeycomb that lights up where you draw
     ===================================================================== */
  (function field() {
    var cv = document.getElementById("hexfield");
    if (!cv || !cv.getContext) return;
    var ctx = cv.getContext("2d"), hexes = [], W = 0, H = 0, R = 14, line, accent, running = false, lastT = 0;
    function colours() { line = tokenRGB("--line", "#d2d2d7"); accent = tokenRGB("--accent", "#e85d04"); }
    function layout() {
      var rect = cv.getBoundingClientRect(), dpr = Math.min(window.devicePixelRatio || 1, 2);
      W = rect.width; H = rect.height; cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      R = W < 500 ? 12 : 15;
      var old = hexes; hexes = [];
      var dx = SQ3 * R, dy = 1.5 * R;
      for (var row = 0, y = R; y < H - R * 0.2; row++, y += dy) {
        for (var x = R + (row % 2 ? dx / 2 : 0); x < W - R * 0.5; x += dx) hexes.push({ x: x, y: y, heat: 0 });
      }
      // A few warm cells so the field doesn't start empty
      for (var k = 0; k < hexes.length; k += 7) hexes[k].base = (k * 37 % 11 === 0) ? 0.55 : 0;
      draw();
    }
    function hexPath(x, y, r) {
      ctx.beginPath();
      for (var v = 0; v < 6; v++) { var a = -PI / 2 + v * PI / 3, X = x + r * Math.cos(a), Y = y + r * Math.sin(a); if (v) ctx.lineTo(X, Y); else ctx.moveTo(X, Y); }
      ctx.closePath();
    }
    function draw() {
      ctx.clearRect(0, 0, W, H);
      var hot = false;
      hexes.forEach(function (h) {
        var heat = Math.max(h.heat, h.base || 0);
        hexPath(h.x, h.y, R * (0.86 + 0.08 * h.heat));
        ctx.fillStyle = css(mixc(line, accent, heat), 0.35 + 0.65 * heat);
        ctx.fill();
        if (h.heat > 0.01) hot = true;
      });
      return hot;
    }
    function loop(now) {
      var dt = Math.min(48, now - (lastT || now)); lastT = now;
      var decay = Math.pow(0.965, dt / 16.7);
      hexes.forEach(function (h) { if (h.delay > 0) { h.delay -= dt; if (h.delay <= 0) h.heat = 1; } else h.heat *= decay; });
      var hot = draw() || hexes.some(function (h) { return h.delay > 0; });
      if (hot) requestAnimationFrame(loop); else { running = false; lastT = 0; }
    }
    function kick() { if (!running) { running = true; requestAnimationFrame(loop); } }
    function local(e) { var r = cv.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; }
    cv.addEventListener("pointermove", function (e) {
      var p = local(e);
      hexes.forEach(function (h) { var d = Math.hypot(h.x - p[0], h.y - p[1]); if (d < R * 2.2) h.heat = Math.max(h.heat, 1 - d / (R * 2.6)); });
      kick();
    }, { passive: true });
    cv.addEventListener("pointerdown", function (e) {   // a tap sends a ripple out across the comb
      var p = local(e);
      hexes.forEach(function (h) { var d = Math.hypot(h.x - p[0], h.y - p[1]); if (d < R * 9) h.delay = reduce ? 1 : d * 3.2; });
      kick();
    }, { passive: true });
    colours(); layout();
    var rt; window.addEventListener("resize", function () { clearTimeout(rt); rt = setTimeout(layout, 150); });
    var mq = window.matchMedia ? window.matchMedia("(prefers-color-scheme: dark)") : null;
    if (mq && mq.addEventListener) mq.addEventListener("change", function () { colours(); draw(); });
  })();

  /* =====================================================================
     4. Therapeutic window: the two bars grow in once, when scrolled into view
     ===================================================================== */
  (function windowChart() {
    var fig = document.getElementById("window");
    if (!fig || reduce || !("IntersectionObserver" in window)) return;
    var r = fig.getBoundingClientRect();
    if (r.top < window.innerHeight && r.bottom > 0) return;   // already on screen: leave it drawn
    fig.classList.add("pre");
    var io = new IntersectionObserver(function (es) {
      es.forEach(function (en) { if (en.isIntersecting) { setTimeout(function () { fig.classList.remove("pre"); }, 80); io.disconnect(); } });
    }, { threshold: 0.45 });
    io.observe(fig);
  })();

  /* =====================================================================
     5. Copy the email address
     ===================================================================== */
  Array.prototype.forEach.call(document.querySelectorAll("[data-copy]"), function (btn) {
    btn.addEventListener("click", function () {
      var text = btn.getAttribute("data-copy");
      function done(label) { btn.textContent = label; setTimeout(function () { btn.textContent = "Copy"; }, 1600); }
      function fallback() {
        var el = document.getElementById("email-text"), range = document.createRange();
        range.selectNodeContents(el); var sel = window.getSelection(); sel.removeAllRanges(); sel.addRange(range);
        done("Selected");
      }
      if (navigator.clipboard && navigator.clipboard.writeText) navigator.clipboard.writeText(text).then(function () { done("Copied"); }, fallback);
      else fallback();
    });
  });
})();
