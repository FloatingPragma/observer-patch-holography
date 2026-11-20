"use strict";
// Each experiment and each OpenWave row has its own page. site.js lists the pages and the
// OpenWave rows; a page states what it shows (window.OPH_PAGE) and embeds its own run.
const SITE = window.OPH_SITE || { runs: [], openwave: {} };
const PAGE = window.OPH_PAGE || { view: "run", root: "./" };
const EMBEDDED = window.OPH_RUNS || [];
const OPENWAVE = SITE.openwave || {};
// A folder address opened from disk needs its index.html.
const href = path => { const url = PAGE.root + path; return location.protocol === "file:" && url.endsWith("/") ? url + "index.html" : url; };
const $ = id => document.getElementById(id);
const C_SIGNAL = 0.5;
const MONO = () => getComputedStyle(document.body).getPropertyValue("--mono");
const PALETTE = ["#4af1f2", "#ffb04a", "#c38bff", "#57d49a", "#ff6b6b"];
const state = { runs: [], run: null, frame: "observer", source: "law", playing: true, t: 0,
  yaw: 0.7, pitch: 0.42, dist: 2.25, last: performance.now() };

// ---------- decoding ----------
async function inflate(bytes, format) {
  const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream(format));
  return new Uint8Array(await new Response(stream).arrayBuffer());
}
const b64 = s => Uint8Array.from(atob(s), c => c.charCodeAt(0));
async function unpack(block) {
  const raw = await inflate(b64(block.zlib_base64), "deflate");
  const Types = { int8: Int8Array, uint8: Uint8Array, float32: Float32Array, float64: Float64Array };
  const T = Types[block.dtype];
  return new T(raw.buffer, raw.byteOffset, raw.byteLength / T.BYTES_PER_ELEMENT);
}
async function openRun(bytes, label) {
  let text;
  try { text = new TextDecoder().decode(await inflate(bytes, "gzip")); }
  catch (e) { text = new TextDecoder().decode(bytes); }
  const doc = JSON.parse(text);
  if (doc.format !== "oph-simulator-run/1") throw new Error("unknown run format " + doc.format);
  await VIEWS[doc.kind].prepare(doc);
  doc.label = (doc.step ? doc.step + " · " : "") + doc.title + (label ? ` (${label})` : "");
  state.runs.push(doc);
  if (label) { doc.local = state.runs.length - 1; $("runPicker").add(new Option(doc.label, "local:" + doc.local)); }
  return doc;
}
function esc(s) { return String(s).replace(/[&<>]/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;" })[c]); }
function fit(cv) {
  const dpr = Math.min(2, window.devicePixelRatio || 1), w = cv.clientWidth, h = cv.clientHeight;
  if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) { cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr); }
  return dpr;
}
function provenance(doc) {
  const p = doc.provenance;
  return `<span>computed</span><b>${esc(p.created)} · ${p.seconds} s</b><span>engine</span><b>Python ${esc(p.python)} · NumPy ${esc(p.numpy)} · ${esc(p.machine)}</b>`;
}

// ---------- WebGL for the field ----------
const canvas = $("gl");
const gl = canvas.getContext("webgl", { antialias: true, premultipliedAlpha: false });
function shader(type, src) {
  const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
  return s;
}
function program(vs, fs) {
  const p = gl.createProgram();
  gl.attachShader(p, shader(gl.VERTEX_SHADER, vs)); gl.attachShader(p, shader(gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(p));
  return p;
}
const pointProg = program(`
attribute vec3 aPos; attribute float aVal;
uniform mat4 uMVP; uniform float uSize; uniform float uPixel;
varying float vVal;
void main() {
  vec4 p = uMVP * vec4(aPos, 1.0);
  gl_Position = p; vVal = aVal;
  float a = aVal > 5.0 ? 0.0 : min(abs(aVal), 1.0);
  float s = aVal > 5.0 ? 0.55 * uSize : uSize * (0.3 + 1.5 * a);
  gl_PointSize = max(1.0, s * uPixel / p.w);
}`, `
precision mediump float;
varying float vVal; uniform float uGain;
void main() {
  vec2 d = gl_PointCoord - 0.5; float r = length(d);
  if (r > 0.5) discard;
  float fall = smoothstep(0.5, 0.05, r);
  if (vVal > 5.0) { gl_FragColor = vec4(vec3(0.23, 0.28, 0.33) * 0.16 * fall, 1.0); return; }
  float a = min(abs(vVal), 1.0);
  vec3 base = vVal >= 0.0 ? vec3(1.0, 0.62, 0.20) : vec3(0.18, 0.58, 1.0);
  vec3 col = mix(base, vec3(1.0), 0.45 * a * a);
  gl_FragColor = vec4(col * pow(a, 1.7) * fall * uGain * 1.5, 1.0);
}`);
const lineProg = program(`
attribute vec3 aPos; uniform mat4 uMVP;
void main() { gl_Position = uMVP * vec4(aPos, 1.0); }`, `
precision mediump float; uniform vec4 uColor;
void main() { gl_FragColor = uColor; }`);
const posBuf = gl.createBuffer(), valBuf = gl.createBuffer(), lineBuf = gl.createBuffer();
function mat4mul(a, b) {
  const o = new Float32Array(16);
  for (let c = 0; c < 4; c++) for (let r = 0; r < 4; r++) { let s = 0; for (let k = 0; k < 4; k++) s += a[k * 4 + r] * b[c * 4 + k]; o[c * 4 + r] = s; }
  return o;
}
function perspective(fov, aspect, near, far) {
  const f = 1 / Math.tan(fov / 2), nf = 1 / (near - far);
  return new Float32Array([f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0, (far + near) * nf, -1, 0, 0, 2 * far * near * nf, 0]);
}
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm3 = a => { const l = Math.hypot(...a); return [a[0] / l, a[1] / l, a[2] / l]; };
function lookAt(eye) {
  const f = norm3([-eye[0], -eye[1], -eye[2]]), s = norm3(cross(f, [0, 1, 0])), u = cross(s, f);
  return new Float32Array([s[0], u[0], -f[0], 0, s[1], u[1], -f[1], 0, s[2], u[2], -f[2], 0, -dot(s, eye), -dot(u, eye), dot(f, eye), 1]);
}

// ---------- view: the 3D field ----------
const UNKNOWN = 9;
const field = {
  footer: true,
  async prepare(doc) {
    const N = doc.carrier.size;
    doc.N = N; doc.N3 = N * N * N; doc.h = 1 / N; doc.F = doc.frames.shape[0];
    doc.duration = doc.chart_time; doc.rate = 0.32;
    doc.values = await unpack(doc.frames);
    doc.pos = new Float32Array(doc.N3 * 3); doc.r = new Float32Array(doc.N3);
    doc.spatial = doc.preparation.modes.map(() => new Float32Array(doc.N3));
    doc.omega = doc.preparation.modes.map(([a, k]) => 0.5 * Math.hypot(2 * Math.PI * Math.hypot(...k), doc.continuum.Q));
    let i = 0, norm = 0;
    for (let x = 0; x < N; x++) for (let y = 0; y < N; y++) for (let z = 0; z < N; z++, i++) {
      const cx = x / N, cy = y / N, cz = z / N;
      doc.pos.set([cx - 0.5, cz - 0.5, cy - 0.5], 3 * i);
      doc.r[i] = Math.hypot(cx - 0.5, cy - 0.5, cz - 0.5);
      let v = 0;
      doc.preparation.modes.forEach(([a, k], j) => { const s = Math.cos(2 * Math.PI * (k[0] * cx + k[1] * cy + k[2] * cz)); doc.spatial[j][i] = s; v += a * s; });
      norm = Math.max(norm, Math.abs(v));
    }
    doc.amps = doc.preparation.modes.map(([a]) => a / norm);
    doc.current = new Float32Array(doc.N3);
  },
  enter(doc) {
    $("legend").innerHTML = `<span><i style="background:var(--warm)"></i>phi &gt; 0</span><span><i style="background:var(--cool)"></i>phi &lt; 0</span>
      <span><i style="background:#3a4652"></i>no record received</span><span><i style="background:var(--accent)"></i>observer A and its record horizon</span>`;
    const rows = doc.convergence.map((c, k) => {
      const prev = doc.convergence[k - 1];
      const order = prev ? Math.log(prev.max_error / c.max_error) / Math.log(c.size / prev.size) : null;
      return `<tr><td>N = ${c.size}</td><td>${c.max_error.toFixed(5)}</td><td>${order ? order.toFixed(2) : ""}</td></tr>`;
    }).join("");
    $("convergence").innerHTML = `<table><tr><th>Carrier</th><th>largest |law − continuum|</th><th>order</th></tr>${rows}</table>
      <p class="note">Deviation of the law's field from the Klein–Gordon solution, relative to the largest initial value, over chart times 0 to ${doc.chart_time}. Order 2 means the error falls by four when the resolution doubles.</p>`;
    $("explain").innerHTML = `The law runs on a carrier of ${doc.N}<sup>3</sup> addresses with a private register in sector Q = 2π·${doc.carrier.private_order}. A field pattern is released at rest. Each point is the scalar field the observers' interferometers read at one address. Observer A sits at the centre. A signal travels at c = 1/2 chart units per chart time, so the picture A holds at its clock time t shows each address as it was at t − r/c. Points that no record has reached are grey. The chart description shows every address at the same chart time: the effective description all observers agree on once every record is collected.`;
    $("meta").innerHTML = `<span>equation</span><b>${esc(doc.continuum.equation)}</b><span>signal speed c</span><b>1/2</b>
      <span>rest frequency</span><b>|Q|/2 = ${doc.continuum.rest_angular_frequency.toFixed(4)}</b><span>updates</span><b>${doc.updates} (h = 1/${doc.N})</b>` + provenance(doc);
  },
  value(doc, i, t) {
    const emitted = state.frame === "observer" ? t - doc.r[i] / C_SIGNAL : t;
    if (emitted < -1e-9) return UNKNOWN;
    const law = () => doc.values[Math.min(doc.F - 1, Math.max(0, Math.round(emitted / doc.h))) * doc.N3 + i] / 127;
    const cont = () => { let v = 0; for (let j = 0; j < doc.amps.length; j++) v += doc.amps[j] * Math.cos(doc.omega[j] * emitted) * doc.spatial[j][i]; return v; };
    if (state.source === "law") return law();
    if (state.source === "continuum") return cont();
    return 20 * (law() - cont());
  },
  lines(doc, t) {
    const L = [], e = 0.5;
    const c = [[-e,-e,-e],[e,-e,-e],[e,e,-e],[-e,e,-e],[-e,-e,e],[e,-e,e],[e,e,e],[-e,e,e]];
    [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]].forEach(([a, b]) => L.push(...c[a], ...c[b]));
    const box = L.length / 3, R = C_SIGNAL * t, M = 96;
    if (state.frame === "observer" && R < 0.88) {
      for (const plane of [[0, 1], [0, 2], [1, 2]]) for (let k = 0; k < M; k++) for (const kk of [k, k + 1]) {
        const p = [0, 0, 0], a = 2 * Math.PI * kk / M; p[plane[0]] = R * Math.cos(a); p[plane[1]] = R * Math.sin(a); L.push(...p);
      }
    }
    const horizon = L.length / 3 - box;
    for (const s of doc.stations) {
      const p = [s[0] / doc.N - 0.5, s[2] / doc.N - 0.5, s[1] / doc.N - 0.5], d = 0.018;
      for (let ax = 0; ax < 3; ax++) { const a = p.slice(), b = p.slice(); a[ax] -= d; b[ax] += d; L.push(...a, ...b); }
    }
    return { data: new Float32Array(L), box, horizon, stations: L.length / 3 - box - horizon };
  },
  draw(doc) {
    fit(canvas);
    gl.viewport(0, 0, canvas.width, canvas.height);
    gl.clearColor(0, 0, 0, 1); gl.clear(gl.COLOR_BUFFER_BIT);
    const eye = [state.dist * Math.cos(state.pitch) * Math.sin(state.yaw), state.dist * Math.sin(state.pitch), state.dist * Math.cos(state.pitch) * Math.cos(state.yaw)];
    const mvp = mat4mul(perspective(0.78, canvas.width / canvas.height, 0.05, 20), lookAt(eye));
    for (let i = 0; i < doc.N3; i++) doc.current[i] = this.value(doc, i, state.t);
    gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ONE); gl.disable(gl.DEPTH_TEST);
    gl.useProgram(pointProg);
    gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
    if (posBuf.owner !== doc) { gl.bufferData(gl.ARRAY_BUFFER, doc.pos, gl.STATIC_DRAW); posBuf.owner = doc; }
    let loc = gl.getAttribLocation(pointProg, "aPos"); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
    gl.bindBuffer(gl.ARRAY_BUFFER, valBuf); gl.bufferData(gl.ARRAY_BUFFER, doc.current, gl.DYNAMIC_DRAW);
    const vloc = gl.getAttribLocation(pointProg, "aVal"); gl.enableVertexAttribArray(vloc); gl.vertexAttribPointer(vloc, 1, gl.FLOAT, false, 0, 0);
    gl.uniformMatrix4fv(gl.getUniformLocation(pointProg, "uMVP"), false, mvp);
    gl.uniform1f(gl.getUniformLocation(pointProg, "uSize"), +$("size").value * 32 / doc.N);
    gl.uniform1f(gl.getUniformLocation(pointProg, "uPixel"), canvas.height / 520);
    gl.uniform1f(gl.getUniformLocation(pointProg, "uGain"), +$("gain").value);
    gl.drawArrays(gl.POINTS, 0, doc.N3);
    gl.disableVertexAttribArray(vloc);
    const lines = this.lines(doc, state.t);
    gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
    gl.useProgram(lineProg);
    gl.bindBuffer(gl.ARRAY_BUFFER, lineBuf); gl.bufferData(gl.ARRAY_BUFFER, lines.data, gl.DYNAMIC_DRAW);
    loc = gl.getAttribLocation(lineProg, "aPos"); gl.enableVertexAttribArray(loc); gl.vertexAttribPointer(loc, 3, gl.FLOAT, false, 0, 0);
    gl.uniformMatrix4fv(gl.getUniformLocation(lineProg, "uMVP"), false, mvp);
    const color = gl.getUniformLocation(lineProg, "uColor");
    gl.uniform4f(color, 0.45, 0.58, 0.66, 0.35); gl.drawArrays(gl.LINES, 0, lines.box);
    gl.uniform4f(color, 0.29, 0.95, 0.95, 0.75); gl.drawArrays(gl.LINES, lines.box, lines.horizon);
    gl.uniform4f(color, 0.29, 0.95, 0.95, 0.95); gl.drawArrays(gl.LINES, lines.box + lines.horizon, lines.stations);
    this.slice(doc); this.series(doc);
    const n = Math.round(state.t / doc.h), R = C_SIGNAL * state.t;
    $("clock").innerHTML = `A's clock t = <b>${state.t.toFixed(3)}</b> <small>chart units · ${n} of ${doc.updates} updates</small>`;
    $("hud").innerHTML = `<b>${state.frame === "observer" ? "Observer A's records" : "Chart description (audit)"}</b> · ${({ law: "the law", continuum: "the continuum limit", difference: "law − continuum, ×20" })[state.source]}<br>` +
      (state.frame === "observer" ? (R < 0.87 ? `record horizon r = c·t = ${R.toFixed(3)}` : "every address has reached A") : "all addresses at one chart time");
  },
  slice(doc) {
    const cv = $("slice"), ctx = cv.getContext("2d"), fine = state.source === "continuum";
    const M = fine ? 140 : doc.N, img = ctx.createImageData(M, M), zc = Math.floor(doc.N / 2);
    for (let a = 0; a < M; a++) for (let b = 0; b < M; b++) {
      let v;
      if (fine) {
        const cx = a / M, cy = b / M, r = Math.hypot(cx - 0.5, cy - 0.5);
        const emitted = state.frame === "observer" ? state.t - r / C_SIGNAL : state.t;
        if (emitted < 0) v = UNKNOWN;
        else { v = 0; doc.preparation.modes.forEach(([amp, k], j) => { v += doc.amps[j] * Math.cos(doc.omega[j] * emitted) * Math.cos(2 * Math.PI * (k[0] * cx + k[1] * cy + k[2] * 0.5)); }); }
      } else v = doc.current[(a * doc.N + b) * doc.N + zc];
      const [r, g, bl] = diverging(v), o = 4 * ((M - 1 - b) * M + a);
      img.data[o] = r; img.data[o + 1] = g; img.data[o + 2] = bl; img.data[o + 3] = 255;
    }
    const off = document.createElement("canvas"); off.width = M; off.height = M; off.getContext("2d").putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = fine; ctx.clearRect(0, 0, cv.width, cv.height); ctx.drawImage(off, 0, 0, cv.width, cv.height);
    const R = C_SIGNAL * state.t * cv.width;
    ctx.strokeStyle = "rgba(74,241,242,.9)"; ctx.lineWidth = 2;
    if (state.frame === "observer" && R < cv.width) { ctx.beginPath(); ctx.arc(cv.width / 2, cv.height / 2, R, 0, 2 * Math.PI); ctx.stroke(); }
    ctx.fillStyle = "#4af1f2"; ctx.beginPath(); ctx.arc(cv.width / 2, cv.height / 2, 6, 0, 2 * Math.PI); ctx.fill();
    ctx.font = "22px " + MONO(); ctx.fillText("A", cv.width / 2 + 10, cv.height / 2 - 10);
    $("sliceNote").textContent = fine
      ? "The continuum limit evaluated on a fine grid: a smooth field." + (state.frame === "observer" ? " Inside the circle A holds records; outside it A holds no record." : "")
      : `The law's ${doc.N} × ${doc.N} addresses in this plane, one value per address.` + (state.frame === "observer" ? " Inside the circle A holds records." : "");
  },
  series(doc) {
    const cv = $("series"), ctx = cv.getContext("2d"), W = cv.width, H = cv.height, pad = 34;
    ctx.clearRect(0, 0, W, H);
    const T = doc.chart_time, X = t => pad + (W - pad - 10) * t / T, Y = v => H / 2 - (H / 2 - 16) * v;
    ctx.strokeStyle = "#1f2c37"; ctx.lineWidth = 1; ctx.beginPath(); ctx.moveTo(pad, Y(0)); ctx.lineTo(W - 10, Y(0)); ctx.stroke();
    ctx.fillStyle = "#8397a6"; ctx.font = "18px " + MONO();
    ctx.fillText("+1", 2, Y(1) + 6); ctx.fillText("−1", 2, Y(-1)); ctx.fillText("chart time of emission →", pad, H - 4);
    for (const [j, col] of [[0, "#4af1f2"], [1, "#ffb04a"], [8, "#c38bff"]]) {
      const s = doc.stations[j], r = Math.hypot(s[0] / doc.N - 0.5, s[1] / doc.N - 0.5, s[2] / doc.N - 0.5);
      const limit = state.frame === "observer" ? state.t - r / C_SIGNAL : state.t;
      ctx.setLineDash([6, 6]); ctx.strokeStyle = col + "66"; ctx.lineWidth = 1.5; ctx.beginPath();
      doc.station_continuum.forEach((row, n) => { const t = n * doc.h; n ? ctx.lineTo(X(t), Y(row[j])) : ctx.moveTo(X(t), Y(row[j])); });
      ctx.stroke(); ctx.setLineDash([]);
      if (limit < 0) continue;
      ctx.strokeStyle = col; ctx.lineWidth = 2.5; ctx.beginPath(); let started = false;
      doc.station_law.forEach((row, n) => { const t = n * doc.h; if (t > limit + 1e-9) return; started ? ctx.lineTo(X(t), Y(row[j])) : ctx.moveTo(X(t), Y(row[j])); started = true; });
      ctx.stroke();
    }
    ctx.fillStyle = "#4af1f2"; ctx.fillText("A", W - 140, 22); ctx.fillStyle = "#ffb04a"; ctx.fillText("corner 1", W - 116, 22); ctx.fillStyle = "#c38bff"; ctx.fillText("8", W - 20, 22);
  },
};
function diverging(v) {
  if (v === UNKNOWN) return [16, 21, 27];
  const a = Math.min(1, Math.abs(v)), base = v >= 0 ? [255, 158, 52] : [46, 148, 255];
  const k = Math.pow(a, 0.85), w = 0.45 * a * a;
  return base.map(c => Math.round((c * (1 - w) + 255 * w) * k + 7 * (1 - k)));
}

// ---------- 2D helpers ----------
const c2d = $("c2d");
function axes(ctx, box, xr, yr, xl, yl, opts = {}) {
  const [x0, y0, w, h] = box, X = v => x0 + w * (v - xr[0]) / (xr[1] - xr[0]), Y = v => y0 + h - h * (v - yr[0]) / (yr[1] - yr[0]);
  ctx.strokeStyle = "#1f2c37"; ctx.lineWidth = 1; ctx.fillStyle = "#8397a6"; ctx.font = `${opts.font || 12}px ${MONO()}`;
  for (const v of opts.yticks || []) { ctx.beginPath(); ctx.moveTo(x0, Y(v)); ctx.lineTo(x0 + w, Y(v)); ctx.stroke(); ctx.textAlign = "right"; ctx.fillText(opts.yfmt ? opts.yfmt(v) : v, x0 - 6, Y(v) + 4); }
  for (const v of opts.xticks || []) { ctx.textAlign = "center"; ctx.fillText(opts.xfmt ? opts.xfmt(v) : v, X(v), y0 + h + 16); }
  ctx.strokeStyle = "#33444f"; ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x0, y0 + h); ctx.lineTo(x0 + w, y0 + h); ctx.stroke();
  ctx.textAlign = "center"; ctx.fillText(xl, x0 + w / 2, y0 + h + 34);
  ctx.save(); ctx.translate(x0 - 46, y0 + h / 2); ctx.rotate(-Math.PI / 2); ctx.fillText(yl, 0, 0); ctx.restore();
  ctx.textAlign = "left";
  return { X, Y };
}
function ticks(lo, hi, n) { const step = Math.pow(10, Math.floor(Math.log10((hi - lo) / n))); const s = [1, 2, 5, 10].map(m => m * step).find(m => (hi - lo) / m <= n); const out = []; for (let v = Math.ceil(lo / s) * s; v <= hi + 1e-12; v += s) out.push(+v.toFixed(10)); return out; }

// ---------- view: a journey with clocks ----------
const NAMES = { beat: "Beat clock (two private modes)", optical: "Optical clock (pulse on a transverse ring)", rails: "Rail clock (two rails of one constituent)" };
const journey = {
  footer: true,
  async prepare(doc) {
    doc.values = await unpack(doc.frames);
    [doc.F, , doc.M] = doc.frames.shape;
    doc.duration = doc.reunion; doc.rate = doc.reunion / 16;
  },
  enter(doc) {
    $("legend").innerHTML = `<span><i style="background:var(--accent)"></i>traveller</span><span><i style="background:var(--warm)"></i>stationary twin</span><span><i style="background:#fff"></i>home laboratory</span>`;
    $("dyn").innerHTML = doc.clocks.map((c, k) => `<section><h2>${esc(NAMES[c.clock] || c.clock)}</h2>
      <div class="dials"><figure><canvas id="dialR${k}" width="300" height="300" class="panel"></canvas><figcaption>stationary twin</figcaption></figure>
      <figure><canvas id="dialM${k}" width="300" height="300" class="panel"></canvas><figcaption>traveller</figcaption></figure></div>
      <canvas id="adv${k}" width="680" height="300" class="panel" style="margin-top:10px"></canvas>
      <div class="kv" id="read${k}"></div></section>`).join("") +
      `<section><h2>At reunion</h2><table><tr><th>clock</th><th>forecast</th><th>reading</th><th>difference</th></tr>${doc.clocks.map(c =>
        `<tr><td>${esc(c.clock)}</td><td>${doc.forecast.toFixed(5)}</td><td class="good">${c.ratio.toFixed(5)}</td><td>${(c.ratio - doc.forecast >= 0 ? "+" : "") + (c.ratio - doc.forecast).toFixed(5)}</td></tr>`).join("")}</table>
      <p class="note">The forecast is the flat metric, sqrt(1 − v²/c²), with v from the traveller's own envelope and c the archived free-wave speed ${doc.reference_speed.toFixed(6)}. The reading is the traveller's advance at home at reunion divided by its twin's advance over the same ${doc.reunion} home ticks.</p></section>`;
    $("explain").innerHTML = esc(doc.summary) + " The picture of the carrier and the traveller's running dial are audit views computed from the law; the reunion readings use only the home laboratory's records.";
    $("meta").innerHTML = `<span>carrier</span><b>${doc.carrier.size} × ${doc.carrier.size} addresses in x, y</b><span>route</span><b>(${doc.route.join(", ")})</b>
      <span>speed</span><b>${doc.speed.toFixed(5)} addresses per tick</b><span>reunion</span><b>${doc.reunion} home ticks</b>` + provenance(doc);
  },
  frameAt(doc) { const ticks = doc.ticks; let f = 0; while (f < ticks.length - 1 && ticks[f + 1] <= state.t) f++; return f; },
  draw(doc) {
    const dpr = fit(c2d), ctx = c2d.getContext("2d"), W = c2d.width, H = c2d.height;
    ctx.fillStyle = "#000"; ctx.fillRect(0, 0, W, H);
    const f = this.frameAt(doc), M = doc.M, S = Math.min(W - 40 * dpr, H - 70 * dpr), x0 = (W - S) / 2, y0 = (H - S) / 2 + 10 * dpr;
    const img = ctx.createImageData(M, M), base = f * 2 * M * M;
    for (let i = 0; i < M; i++) for (let j = 0; j < M; j++) {
      const a = doc.values[base + i * M + j] / 255, b = doc.values[base + M * M + i * M + j] / 255;
      const u = (i + M / 2) % M, v = (j + M / 2) % M, o = 4 * ((M - 1 - v) * M + u);
      const c = 0.62 * b;
      img.data[o] = Math.min(255, 7 + 74 * a + 255 * c); img.data[o + 1] = Math.min(255, 10 + 241 * a + 176 * c);
      img.data[o + 2] = Math.min(255, 14 + 242 * a + 74 * c); img.data[o + 3] = 255;
    }
    const off = document.createElement("canvas"); off.width = M; off.height = M; off.getContext("2d").putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true; ctx.drawImage(off, x0, y0, S, S);
    ctx.strokeStyle = "#1f2c37"; ctx.lineWidth = 1 * dpr; ctx.strokeRect(x0, y0, S, S);
    const N = doc.carrier.size, P = ([x, y]) => [x0 + S * (((x + N / 2) % N) / N), y0 + S - S * (((y + N / 2) % N) / N)];
    ctx.fillStyle = "rgba(74,241,242,.55)";
    for (let k = 0; k <= f; k += 2) { const [px, py] = P(doc.centre[k]); ctx.fillRect(px - dpr, py - dpr, 2 * dpr, 2 * dpr); }
    const [hx, hy] = P([0, 0]);
    ctx.strokeStyle = "#ffffff"; ctx.lineWidth = 1.5 * dpr; ctx.beginPath(); ctx.arc(hx, hy, 9 * dpr, 0, 2 * Math.PI); ctx.stroke();
    ctx.fillStyle = "#dce7ee"; ctx.font = `${12 * dpr}px ${MONO()}`; ctx.fillText("home", hx + 12 * dpr, hy - 10 * dpr);
    const tick = doc.ticks[f];
    $("hud").innerHTML = `<b>${esc(doc.step)} · ${esc(doc.title)}</b><br>periodic carrier, home at the centre · audit view of the law's probabilities`;
    $("clock").innerHTML = `home clock <b>${tick}</b> <small>of ${doc.reunion} ticks${tick === doc.reunion ? " · reunion" : ""}</small>`;
    doc.clocks.forEach((c, k) => this.dials(doc, c, k, f));
  },
  dials(doc, c, k, f) {
    const total = c.rest[c.rest.length - 1], sweep = 2 * Math.PI * 0.9;
    const dial = (id, value, color, mark) => {
      const cv = $(id); if (!cv) return; const ctx = cv.getContext("2d"), W = cv.width, R = W * 0.4, cx = W / 2, cy = W / 2;
      ctx.clearRect(0, 0, W, W);
      ctx.strokeStyle = "#1f2c37"; ctx.lineWidth = 10; ctx.beginPath(); ctx.arc(cx, cy, R, -Math.PI / 2, -Math.PI / 2 + sweep); ctx.stroke();
      ctx.strokeStyle = color; ctx.beginPath(); ctx.arc(cx, cy, R, -Math.PI / 2, -Math.PI / 2 + sweep * Math.max(0, value / total)); ctx.stroke();
      if (mark !== undefined) { const a = -Math.PI / 2 + sweep * mark; ctx.strokeStyle = "#ffffff"; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(cx + (R - 16) * Math.cos(a), cy + (R - 16) * Math.sin(a)); ctx.lineTo(cx + (R + 16) * Math.cos(a), cy + (R + 16) * Math.sin(a)); ctx.stroke(); }
      const a = -Math.PI / 2 + sweep * Math.max(0, value / total);
      ctx.fillStyle = "#dce7ee"; ctx.beginPath(); ctx.arc(cx + R * Math.cos(a), cy + R * Math.sin(a), 9, 0, 2 * Math.PI); ctx.fill();
      ctx.fillStyle = "#dce7ee"; ctx.font = `34px ${MONO()}`; ctx.textAlign = "center"; ctx.fillText((value / total).toFixed(3), cx, cy + 12);
      ctx.fillStyle = "#8397a6"; ctx.font = `17px ${MONO()}`; ctx.fillText("of the twin's total", cx, cy + 40); ctx.textAlign = "left";
    };
    dial(`dialR${k}`, c.rest[f], "#ffb04a");
    dial(`dialM${k}`, c.moving[f], "#4af1f2", doc.forecast);
    const cv = $(`adv${k}`), ctx = cv.getContext("2d"), W = cv.width, H = cv.height;
    ctx.clearRect(0, 0, W, H);
    const { X, Y } = axes(ctx, [64, 14, W - 84, H - 64], [0, doc.reunion], [0, 1.05], "home ticks", "advance / twin's total",
      { yticks: [0, 0.5, 1], xticks: ticks(0, doc.reunion, 4), font: 16 });
    ctx.setLineDash([8, 6]); ctx.strokeStyle = "#ffffff88"; ctx.lineWidth = 2; ctx.beginPath(); ctx.moveTo(X(0), Y(0)); ctx.lineTo(X(doc.reunion), Y(doc.forecast)); ctx.stroke(); ctx.setLineDash([]);
    for (const [arr, col] of [[c.rest, "#ffb04a"], [c.moving, "#4af1f2"]]) {
      ctx.strokeStyle = col; ctx.lineWidth = 3; ctx.beginPath();
      for (let i = 0; i <= f; i++) { const x = X(doc.ticks[i]), y = Y(arr[i] / total); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }
      ctx.stroke();
    }
    ctx.strokeStyle = "#4b5d6b"; ctx.beginPath(); ctx.moveTo(X(doc.ticks[f]), Y(0)); ctx.lineTo(X(doc.ticks[f]), Y(1.05)); ctx.stroke();
    const running = c.rest[f] > 0 ? c.moving[f] / c.rest[f] : NaN;
    $(`read${k}`).innerHTML = `<span>twin</span><b>${c.rest[f].toFixed(4)} ${c.unit}</b><span>traveller</span><b>${c.moving[f].toFixed(4)} ${c.unit}</b>
      <span>running ratio</span><b>${isFinite(running) ? running.toFixed(4) : "–"}</b><span>white mark</span><b>forecast ${doc.forecast.toFixed(5)}</b>`;
  },
};

// ---------- view: the signal race ----------
const race = {
  footer: true,
  async prepare(doc) {
    doc.values = await unpack(doc.strip);
    [doc.F, doc.L] = doc.strip.shape;
    doc.duration = doc.deadlines.late; doc.rate = doc.deadlines.late / 14;
  },
  enter(doc) {
    const d = doc.deadlines;
    $("legend").innerHTML = `<span><i style="background:#2c6f74"></i>exponentially small precursor</span><span><i style="background:#fff"></i>the wave</span><span><i style="background:var(--warm)"></i>smooth cone, c = 1/2</span><span><i style="background:var(--bad)"></i>relay, excluded access</span>`;
    $("dyn").innerHTML = `<section><h2>What B's detector records</h2><table>
      <tr><td>first nonzero probability</td><td>update ${doc.first_nonzero}</td></tr>
      <tr><td>largest before the early deadline</td><td>${doc.largest_before_early.toExponential(2)}</td></tr>
      <tr><td>smooth arrival, r / c</td><td>update ${d.smooth}</td></tr>
      <tr><td>summed detection probability by the late deadline</td><td>${doc.received_by_late.toFixed(3)}</td></tr></table>
      <p class="note">The law's strict front moves one address per update, so B's probability is exactly zero before update ${d.strict_front}. Between the strict front and the smooth cone each further forward hop costs a factor 1/16, and the probability stays below ${doc.largest_before_early.toExponential(1)} until the early deadline at update ${d.early}. A field-only sender therefore cannot get a fresh bit to B ahead of the smooth cone. A relay that inhabitants could trigger would reach B at update ${d.strict_front} with certainty; that access rule is excluded.</p></section>`;
    $("explain").innerHTML = esc(doc.summary);
    $("meta").innerHTML = `<span>carrier</span><b>a strip of ${doc.carrier.size} addresses</b><span>A, B</span><b>${doc.source}, ${doc.receiver} (distance ${doc.distance})</b>` + provenance(doc);
  },
  draw(doc) {
    const dpr = fit(c2d), ctx = c2d.getContext("2d"), W = c2d.width, H = c2d.height, d = doc.deadlines;
    ctx.fillStyle = "#000"; ctx.fillRect(0, 0, W, H);
    const f = Math.min(doc.F - 1, Math.floor(state.t / doc.stride)), t = f * doc.stride;
    const lo = doc.source - doc.distance * 1.5, hi = doc.receiver + doc.distance * 0.75;
    const x0 = 60 * dpr, w = W - 100 * dpr, X = x => x0 + w * (x - lo) / (hi - lo);
    const top = 70 * dpr, sh = Math.min(120 * dpr, H * 0.22);
    for (let x = Math.ceil(lo); x < hi; x++) {
      const code = doc.values[f * doc.L + x], logp = code / 255 * 60 - 60;
      const faint = code > 0 ? 0.10 + 0.28 * (logp + 60) / 54 : 0, bright = Math.max(0, Math.min(1, (logp + 6) / 4));
      ctx.fillStyle = `rgb(${Math.round(6 + 40 * faint + 245 * bright)},${Math.round(9 + 190 * faint + 240 * bright)},${Math.round(13 + 200 * faint + 238 * bright)})`;
      ctx.fillRect(X(x), top, Math.ceil(w / (hi - lo)) + 1, sh);
    }
    const mark = (x, color, label, dash) => { ctx.strokeStyle = color; ctx.setLineDash(dash || []); ctx.lineWidth = 2 * dpr; ctx.beginPath(); ctx.moveTo(X(x), top - 8 * dpr); ctx.lineTo(X(x), top + sh + 8 * dpr); ctx.stroke(); ctx.setLineDash([]); if (label) { ctx.fillStyle = color; ctx.font = `${12 * dpr}px ${MONO()}`; ctx.fillText(label, X(x) - 6 * dpr, top - 14 * dpr); } };
    mark(doc.source, "#ffffff", "A"); mark(doc.receiver, "#ffffff", "B");
    mark(doc.source + C_SIGNAL * t, "#ffb04a", null, [6, 4]); mark(doc.source - C_SIGNAL * t, "#ffb04a", null, [6, 4]);
    if (t <= doc.distance) { ctx.strokeStyle = "#ff6b6b"; ctx.lineWidth = 2 * dpr; ctx.beginPath(); ctx.arc(X(doc.source + t), top + sh / 2, 9 * dpr, 0, 2 * Math.PI); ctx.stroke(); }
    ctx.fillStyle = "#8397a6"; ctx.font = `${12 * dpr}px ${MONO()}`; ctx.fillText("address →", x0, top + sh + 26 * dpr);
    const cy = top + sh + 60 * dpr, ch = H - cy - 90 * dpr;
    const { X: TX, Y: TY } = axes(ctx, [x0, cy, w, ch], [0, d.late], [-60, 0], "update at B", "log10 P(detect at B)",
      { yticks: [-60, -45, -30, -15, 0], xticks: [0, d.strict_front, d.early, d.smooth, d.late], font: 12 * dpr });
    for (const [x, col, label] of [[d.strict_front, "#ff6b6b", "strict front"], [d.early, "#8397a6", "early deadline"], [d.smooth, "#ffb04a", "smooth arrival"], [d.late, "#8397a6", "late deadline"]]) {
      ctx.strokeStyle = col; ctx.setLineDash([4, 4]); ctx.lineWidth = 1 * dpr; ctx.beginPath(); ctx.moveTo(TX(x), cy); ctx.lineTo(TX(x), cy + ch); ctx.stroke(); ctx.setLineDash([]);
      ctx.fillStyle = col; ctx.fillText(label, TX(x) + 4 * dpr, cy + 14 * dpr);
    }
    ctx.strokeStyle = "#4af1f2"; ctx.lineWidth = 2.5 * dpr; ctx.beginPath(); let on = false;
    for (let n = 0; n <= Math.min(t, d.late); n++) {
      const p = doc.detector[n]; if (!(p > 0)) { on = false; continue; }
      const y = TY(Math.max(-60, Math.log10(p))); on ? ctx.lineTo(TX(n), y) : ctx.moveTo(TX(n), y); on = true;
    }
    ctx.stroke();
    $("hud").innerHTML = `<b>${esc(doc.step)} · ${esc(doc.title)}</b><br>field-only access: the free wave is the only messenger`;
    $("clock").innerHTML = `A's clock <b>${t}</b> <small>updates since the launch</small>`;
  },
};

// ---------- view: a scan ----------
const scan = {
  footer: false,
  async prepare(doc) { doc.duration = 0; },
  enter(doc) {
    $("legend").innerHTML = "";
    const [head, ...rows] = doc.table;
    $("dyn").innerHTML = `<section><h2>Values</h2><table><tr>${head.map(h => `<th>${esc(h)}</th>`).join("")}</tr>${rows.map(r => `<tr>${r.map(v => `<td>${v === null ? "–" : esc(v)}</td>`).join("")}</tr>`).join("")}</table></section>`;
    $("explain").innerHTML = esc(doc.summary);
    $("meta").innerHTML = provenance(doc);
    $("hud").innerHTML = `<b>${esc(doc.step)} · ${esc(doc.title)}</b>`;
  },
  draw(doc) {
    const dpr = fit(c2d), ctx = c2d.getContext("2d"), W = c2d.width, H = c2d.height;
    ctx.fillStyle = "#000"; ctx.fillRect(0, 0, W, H);
    const all = doc.series.flatMap(s => s.points).concat(doc.curves.flatMap(c => c.points)).filter(p => p[1] !== null);
    const xs = all.map(p => p[0]), ys = all.map(p => p[1]);
    const padx = (Math.max(...xs) - Math.min(...xs)) * 0.06 || 1, pady = (Math.max(...ys) - Math.min(...ys)) * 0.08 || 0.1;
    const xr = [Math.min(...xs) - padx, Math.max(...xs) + padx], yr = [Math.min(...ys) - pady, Math.max(...ys) + pady];
    const box = [90 * dpr, 70 * dpr, W - 130 * dpr, H - 160 * dpr];
    const { X, Y } = axes(ctx, box, xr, yr, doc.x_label, doc.y_label, { xticks: ticks(xr[0], xr[1], 6), yticks: ticks(yr[0], yr[1], 6), font: 13 * dpr, yfmt: v => +v.toFixed(4) });
    if (doc.zero_line) { ctx.strokeStyle = "#dce7ee55"; ctx.lineWidth = 1.5 * dpr; ctx.beginPath(); ctx.moveTo(X(xr[0]), Y(0)); ctx.lineTo(X(xr[1]), Y(0)); ctx.stroke(); }
    doc.curves.forEach(c => { ctx.strokeStyle = "#dce7eeaa"; ctx.lineWidth = 2 * dpr; ctx.setLineDash([8, 6]); ctx.beginPath(); c.points.forEach((p, i) => i ? ctx.lineTo(X(p[0]), Y(p[1])) : ctx.moveTo(X(p[0]), Y(p[1]))); ctx.stroke(); ctx.setLineDash([]); });
    doc.series.forEach((s, k) => {
      const col = PALETTE[k % PALETTE.length];
      ctx.strokeStyle = col + "88"; ctx.lineWidth = 1.5 * dpr; ctx.beginPath(); let on = false;
      s.points.forEach(p => { if (p[1] === null) return; on ? ctx.lineTo(X(p[0]), Y(p[1])) : ctx.moveTo(X(p[0]), Y(p[1])); on = true; }); ctx.stroke();
      s.points.forEach(p => { if (p[1] === null) return; ctx.fillStyle = col; ctx.beginPath(); ctx.arc(X(p[0]), Y(p[1]), (5 - k * 0.6) * dpr, 0, 2 * Math.PI); ctx.fill(); });
    });
    ctx.font = `${13 * dpr}px ${MONO()}`;
    const legend = doc.series.map((s, k) => [s.name, PALETTE[k % PALETTE.length], false]).concat(doc.curves.map(c => [c.name, "#dce7ee", true]));
    legend.forEach(([name, col, dashed], k) => {
      const lx = box[0] + 12 * dpr, ly = box[1] + 18 * dpr + k * 20 * dpr;
      ctx.strokeStyle = col; ctx.fillStyle = col; ctx.lineWidth = 2 * dpr;
      if (dashed) { ctx.setLineDash([6, 4]); ctx.beginPath(); ctx.moveTo(lx, ly - 4 * dpr); ctx.lineTo(lx + 22 * dpr, ly - 4 * dpr); ctx.stroke(); ctx.setLineDash([]); }
      else { ctx.beginPath(); ctx.arc(lx + 11 * dpr, ly - 4 * dpr, 5 * dpr, 0, 2 * Math.PI); ctx.fill(); }
      ctx.fillStyle = "#dce7ee"; ctx.fillText(name, lx + 30 * dpr, ly);
    });
  },
};

// ---------- view: formation through the capture vertex ----------
const formation = {
  footer: true,
  async prepare(doc) {
    doc.values = await unpack(doc.joint);
    [doc.F, doc.N] = doc.joint.shape;
    doc.duration = doc.updates; doc.rate = doc.updates / 12;
  },
  enter(doc) {
    $("legend").innerHTML = `<span><i style="background:var(--accent)"></i>the two separate constituents</span><span><i style="background:var(--warm)"></i>the formed pair</span><span><i style="background:#c38bff"></i>markers on the junction wires</span>`;
    $("dyn").innerHTML = `<section><h2>Formation</h2><table>
      <tr><td>probability of a formed pair</td><td>${doc.captured[doc.captured.length - 1].toFixed(4)}</td></tr>
      <tr><td>with the vertex disabled</td><td>${doc.captured_without_vertex.toFixed(4)}</td></tr>
      <tr><td>formed pair ahead of its birthplace, markers traced out</td><td>${(100 * doc.forward_fraction).toFixed(1)}%</td></tr>
      <tr><td>the same, markers read at a coherent collector</td><td>${(100 * doc.forward_fraction_collected).toFixed(1)}%</td></tr></table>
      <p class="note">Once formed, the pair never separates: the contact law keeps its bound subspace invariant. Each marker records the junction and the update of a birth, so with the markers traced out every birth is a localized pair that spreads both ways. A collector that joins all junction wires coherently leaves the birthplace unresolved at its common port; conditioned on that port the pair keeps the incoming net momentum, and the opening of its hop sends it on in that direction.</p>
      <canvas id="capCurve" width="680" height="280" class="panel"></canvas></section>`;
    $("explain").innerHTML = esc(doc.summary) + " The pictures are audit views of the law's probabilities.";
    $("meta").innerHTML = `<span>carrier</span><b>a strip of ${doc.carrier.size} junctions</b><span>launch orders</span><b>${doc.orders.join(", ")}</b>
      <span>conversion angle</span><b>${doc.angle.toFixed(4)} rad</b><span>contact phase</span><b>${doc.contact_phase.toFixed(6)}</b>` + provenance(doc);
  },
  draw(doc) {
    const dpr = fit(c2d), ctx = c2d.getContext("2d"), W = c2d.width, H = c2d.height, N = doc.N;
    ctx.fillStyle = "#000"; ctx.fillRect(0, 0, W, H);
    const f = Math.min(doc.F - 1, Math.round(state.t)), font = `${12 * dpr}px ${MONO()}`;
    const S = Math.min(W * 0.42, H - 120 * dpr), x0 = 50 * dpr, y0 = 70 * dpr;
    const img = ctx.createImageData(N, N);
    for (let a = 0; a < N; a++) for (let b = 0; b < N; b++) {
      const v = doc.values[f * N * N + a * N + b] / 255, o = 4 * ((N - 1 - b) * N + a);
      img.data[o] = 7 + 67 * v; img.data[o + 1] = 10 + 231 * v; img.data[o + 2] = 14 + 228 * v; img.data[o + 3] = 255;
    }
    const off = document.createElement("canvas"); off.width = N; off.height = N; off.getContext("2d").putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = false; ctx.drawImage(off, x0, y0, S, S);
    ctx.strokeStyle = "#ffb04a88"; ctx.setLineDash([6, 5]); ctx.lineWidth = 1.5 * dpr; ctx.beginPath(); ctx.moveTo(x0, y0 + S); ctx.lineTo(x0 + S, y0); ctx.stroke(); ctx.setLineDash([]);
    ctx.fillStyle = "#8397a6"; ctx.font = font;
    ctx.fillText("position of A →", x0, y0 + S + 18 * dpr);
    ctx.save(); ctx.translate(x0 - 10 * dpr, y0 + S); ctx.rotate(-Math.PI / 2); ctx.fillText("position of B →", 0, 0); ctx.restore();
    ctx.fillText("the two constituents while unjoined (dashed: same junction)", x0, y0 - 14 * dpr);
    const sx = x0 + S + 70 * dpr, sw = W - sx - 30 * dpr, base = y0 + S, Z = z => sx + sw * (z + 0.5) / N;
    const wireTop = y0 + 20 * dpr, wireH = base - wireTop - 120 * dpr;
    ctx.strokeStyle = "#1f2c37"; ctx.lineWidth = 1 * dpr;
    for (let z = 0; z < N; z++) { ctx.beginPath(); ctx.moveTo(Z(z), wireTop); ctx.lineTo(Z(z), wireTop + wireH); ctx.stroke(); }
    const ages = Math.max(1, doc.updates);
    for (const [z, n, w] of doc.markers) {
      if (n >= f) continue;
      const age = f - n, y = wireTop + wireH - wireH * age / ages, a = Math.min(1, Math.sqrt(w) * 12);
      ctx.fillStyle = `rgba(195,139,255,${a.toFixed(3)})`; ctx.fillRect(Z(z) - 2 * dpr, y - 1.5 * dpr, 4 * dpr, 3 * dpr);
    }
    ctx.fillStyle = "#8397a6"; ctx.fillText("marker wires: a marker climbs one link per update", sx, wireTop - 8 * dpr);
    const plotTop = wireTop + wireH + 20 * dpr, plotH = base - plotTop;
    const bound = doc.bound[f], coll = doc.collected[f];
    const peak = Math.max(1e-12, ...doc.bound[doc.F - 1], ...doc.bound[f]);
    const cpeak = Math.max(1e-12, ...doc.collected[doc.F - 1]);
    const area = (arr, top, col, fill) => { ctx.beginPath(); ctx.moveTo(Z(0), plotTop + plotH); arr.forEach((v, z) => ctx.lineTo(Z(z), plotTop + plotH - plotH * v / top)); ctx.lineTo(Z(N - 1), plotTop + plotH); ctx.closePath(); ctx.fillStyle = fill; ctx.fill(); ctx.strokeStyle = col; ctx.lineWidth = 2 * dpr; ctx.stroke(); };
    area(bound, peak, "#ffb04a", "rgba(255,176,74,.25)");
    ctx.setLineDash([5, 4]); ctx.strokeStyle = "#4af1f2"; ctx.lineWidth = 2 * dpr; ctx.beginPath();
    coll.forEach((v, z) => { const y = plotTop + plotH - plotH * v / cpeak; z ? ctx.lineTo(Z(z), y) : ctx.moveTo(Z(z), y); }); ctx.stroke(); ctx.setLineDash([]);
    ctx.strokeStyle = "#ffffff55"; ctx.beginPath(); ctx.moveTo(Z(doc.birth_centre), plotTop); ctx.lineTo(Z(doc.birth_centre), plotTop + plotH); ctx.stroke();
    ctx.fillStyle = "#8397a6"; ctx.fillText("formed pair along the strip (dashed: given the collector's common port)", sx, plotTop - 6 * dpr);
    ctx.fillText("junction →", sx, base + 18 * dpr);
    const cv = $("capCurve");
    if (cv) {
      const c = cv.getContext("2d"); c.clearRect(0, 0, cv.width, cv.height);
      const top = Math.max(0.05, doc.captured[doc.captured.length - 1] * 1.15);
      const { X, Y } = axes(c, [70, 14, cv.width - 90, cv.height - 64], [0, doc.updates], [0, top], "update", "P(formed pair)", { yticks: ticks(0, top, 4), xticks: ticks(0, doc.updates, 4), font: 16, yfmt: v => +v.toFixed(3) });
      c.strokeStyle = "#ff6b6b"; c.setLineDash([6, 5]); c.lineWidth = 2; c.beginPath(); c.moveTo(X(0), Y(doc.captured_without_vertex)); c.lineTo(X(doc.updates), Y(doc.captured_without_vertex)); c.stroke(); c.setLineDash([]);
      c.strokeStyle = "#ffb04a"; c.lineWidth = 3; c.beginPath(); for (let n = 0; n <= f; n++) n ? c.lineTo(X(n), Y(doc.captured[n])) : c.moveTo(X(n), Y(doc.captured[n])); c.stroke();
      c.fillStyle = "#ff6b6b"; c.font = `15px ${MONO()}`; c.fillText("vertex disabled", X(doc.updates * 0.55), Y(0) - 8);
    }
    $("hud").innerHTML = `<b>${esc(doc.step)} · ${esc(doc.title)}</b><br>capture vertex with an outgoing marker · audit view of the law's probabilities`;
    $("clock").innerHTML = `update <b>${f}</b> <small>of ${doc.updates}</small>`;
  },
};

// ---------- view: a matter-wave relation (OpenWave Row 30) ----------
function hsv(h, s, v) {
  const i = Math.floor(h * 6), f = h * 6 - i, p = v * (1 - s), q = v * (1 - f * s), u = v * (1 - (1 - f) * s);
  return [[v, u, p], [q, v, p], [p, v, u], [p, q, v], [u, p, v], [v, p, q]][((i % 6) + 6) % 6];
}
const dispersion = {
  footer: true,
  async prepare(doc) {
    doc.values = await unpack(doc.frames);
    [doc.F, , , doc.M] = doc.frames.shape;
    doc.duration = doc.updates; doc.rate = doc.updates / 10;
  },
  enter(doc) {
    $("legend").innerHTML = `<span><i style="background:conic-gradient(red,yellow,lime,cyan,blue,magenta,red)"></i>phase of the matter wave</span><span><i style="background:var(--accent)"></i>OPH law</span><span><i style="background:var(--bad)"></i>changed law</span>`;
    const row = rowOf(doc.openwave_row);
    const per = [64, 128, 256].map(P => {
      const c = doc.rows.filter(r => r.law === "oph" && r.period === P), x = doc.rows.filter(r => r.law === "changed" && r.period === P);
      const rng = (a, k, d) => `${Math.min(...a.map(r => r[k])).toFixed(d)}–${Math.max(...a.map(r => r[k])).toFixed(d)}`;
      return `<tr><td>${P}</td><td>${c[0].k.toFixed(4)}</td><td>${rng(c, "kinetic", 4)}</td><td class="good">${c.filter(r => r.passes).length}/14</td><td>${x.filter(r => r.passes).length}/14</td></tr>`;
    }).join("");
    const rec = doc.recorded;
    $("dyn").innerHTML = `<section><h2>OpenWave Row 30</h2>
      <p class="note" style="margin-top:0">OpenWave's test: <b style="color:var(--ink)">${esc(row ? row.test : "emergent ω² = k² + m² for matter waves")}</b>.</p>
      ${rowLinks(row)}
      <table><tr><th>period</th><th>|k|</th><th>kinetic ratio</th><th>OPH law</th><th>changed law</th></tr>${per}</table>
      <p class="note">ω from the native step of the formed pair, fourteen directions per period. Kinetic ratio (ω² − m²)/(c²|k|²) with m = ${doc.rest.toFixed(6)} and c² = ${doc.c2}; the test requires it within ±${(100 * doc.tolerance).toFixed(1)}% of one. The long-wave expansion gives ω² = m² + 2w m cot(m) |k|², coefficient ${doc.coefficient.toFixed(8)}; the closed form agrees with the native step to ${doc.identity_error.toExponential(1)}.</p></section>
      <section><h2>Recorded in the Row 30 experiment</h2><table>${(rec.items || []).map(([k, v]) => `<tr><td>${esc(k)}</td><td>${esc(v)}</td></tr>`).join("")}</table>
      <p class="note">The observers measured rest lines, wave vectors and frequencies from their own records; the counts are sampled from the law. </p></section>`;
    $("explain").innerHTML = esc(doc.summary);
    $("meta").innerHTML = `<span>contact phase</span><b>${doc.phi.toFixed(6)} (m = φ/2)</b>` + provenance(doc);
  },
  draw(doc) {
    const dpr = fit(c2d), ctx = c2d.getContext("2d"), W = c2d.width, H = c2d.height, M = doc.M;
    ctx.fillStyle = "#000"; ctx.fillRect(0, 0, W, H);
    const f = Math.min(doc.F - 1, Math.round(state.t / doc.stride)), font = `${12 * dpr}px ${MONO()}`;
    const gap = 30 * dpr, P = Math.min((W - 3 * gap) / 2, H * 0.42), y0 = 64 * dpr;
    const labels = ["at rest: the phase turns at the rest frequency m", "moving along x: phase fronts travel, the envelope follows"];
    for (let k = 0; k < 2; k++) {
      const img = ctx.createImageData(M, M), base = (f * 2 + k) * 2 * M * M;
      for (let x = 0; x < M; x++) for (let y = 0; y < M; y++) {
        const a = doc.values[base + x * M + y] / 255, ph = doc.values[base + M * M + x * M + y] / 255;
        const [r, g, b] = hsv(ph, 0.85, Math.min(1, a * 1.15)), o = 4 * ((M - 1 - y) * M + x);
        img.data[o] = 255 * r; img.data[o + 1] = 255 * g; img.data[o + 2] = 255 * b; img.data[o + 3] = 255;
      }
      const off = document.createElement("canvas"); off.width = M; off.height = M; off.getContext("2d").putImageData(img, 0, 0);
      const x0 = gap + k * (P + gap);
      ctx.imageSmoothingEnabled = true; ctx.drawImage(off, x0, y0, P, P);
      ctx.strokeStyle = "#1f2c37"; ctx.strokeRect(x0, y0, P, P);
      ctx.fillStyle = "#8397a6"; ctx.font = font; ctx.fillText(labels[k], x0, y0 - 10 * dpr);
    }
    const cy = y0 + P + 50 * dpr, ch = H - cy - 50 * dpr, cw = W - 2 * gap - 60 * dpr;
    const k2max = Math.max(...doc.rows.map(r => r.k * r.k)) * 1.08, m2 = doc.rest * doc.rest;
    const om2max = Math.max(...doc.rows.map(r => r.omega * r.omega)) * 1.06;
    const { X, Y } = axes(ctx, [gap + 60 * dpr, cy, cw, ch], [0, k2max], [m2 * 0.9, om2max], "|k|² (radar units)", "ω²", { xticks: ticks(0, k2max, 5), yticks: ticks(m2 * 0.9, om2max, 4), font: 12 * dpr, yfmt: v => +v.toFixed(4), xfmt: v => +v.toFixed(3) });
    ctx.fillStyle = "rgba(74,241,242,.10)"; ctx.beginPath();
    ctx.moveTo(X(0), Y(m2)); ctx.lineTo(X(k2max), Y(m2 + doc.c2 * k2max * (1 + doc.tolerance))); ctx.lineTo(X(k2max), Y(m2 + doc.c2 * k2max * (1 - doc.tolerance))); ctx.closePath(); ctx.fill();
    ctx.setLineDash([7, 5]); ctx.strokeStyle = "#dce7eecc"; ctx.lineWidth = 1.5 * dpr; ctx.beginPath(); ctx.moveTo(X(0), Y(m2)); ctx.lineTo(X(k2max), Y(m2 + doc.c2 * k2max)); ctx.stroke(); ctx.setLineDash([]);
    for (const r of doc.rows) { ctx.fillStyle = r.law === "oph" ? "#4af1f2" : "#ff6b6b"; ctx.beginPath(); ctx.arc(X(r.k * r.k), Y(r.omega * r.omega), 4 * dpr, 0, 2 * Math.PI); ctx.fill(); }
    ctx.fillStyle = "#dce7ee"; ctx.font = font;
    ctx.fillText("dashed: ω² = m² + c²|k|², band ±3.5% on the kinetic part", X(k2max * 0.02), Y(om2max) + 16 * dpr);
    $("hud").innerHTML = `<b>${esc(doc.step)} · ${esc(doc.title)}</b><br>formed matter of the OPH law · audit view of the law's amplitudes`;
    $("clock").innerHTML = `update <b>${f * doc.stride}</b> <small>of ${doc.updates}</small>`;
  },
};

// ---------- view: persistence of a bound pair (OpenWave Row 4) ----------
const persistence = {
  footer: true,
  async prepare(doc) {
    doc.values = await unpack(doc.frames);
    [doc.F, , doc.N] = doc.frames.shape;
    doc.duration = doc.updates; doc.rate = doc.updates / 12;
  },
  enter(doc) {
    $("legend").innerHTML = `<span><i style="background:var(--good)"></i>the OPH law</span><span><i style="background:var(--bad)"></i>independent transport (control)</span>`;
    const rec = doc.recorded, row = rowOf(4);
    $("dyn").innerHTML = `<section><h2>OpenWave Row 4</h2>
      <p class="note" style="margin-top:0">OpenWave's test: <b style="color:var(--ink)">${esc(row ? row.test : "the particle state persists in free evolution")}</b>.</p>
      ${rowLinks(row)}
      <table><tr><td>eigenvalue of every update, all ages</td><td>−exp(iφ/2)</td></tr>
      <tr><td>largest change of the profile</td><td>${doc.profile_change.toExponential(1)}</td></tr>
      <tr><td>capture projection at every address</td><td>${doc.capture_projection.toExponential(1)}</td></tr>
      <tr><td>matched weight under the law</td><td>${Math.min(...doc.contact_matched).toFixed(6)}</td></tr>
      <tr><td>matched weight, control, after one update</td><td>${doc.control_matched[1].toFixed(6)}</td></tr></table>
      <p class="note">On the complement of the equal-slot vector the contact coin is the number −exp(iφ/2), and the routing swaps the two ends of every edge. The signed circulation lives entirely there, so each update multiplies it by that number, at every age of the pair. A returning marker releases a pair only through its overlap with the equal-slot pair, which is zero at every address.</p></section>
      <section><h2>Recorded in the Row 4 experiment</h2><table>
      <tr><td>cohorts of 4,194,304 trials</td><td>${rec.cohorts}</td></tr><tr><td>separations</td><td>${rec.separations}</td></tr>
      <tr><td>independent transport, separated</td><td>${rec.control_separated[0].toLocaleString("en")} of ${rec.control_separated[1].toLocaleString("en")}</td></tr>
      <tr><td>ordinary pair, first closed return</td><td>${(100 * rec.ordinary_return_separation).toFixed(5)}% separated</td></tr>
      <tr><td>the same with the route cut</td><td>${(100 * rec.route_cut_return_separation).toFixed(5)}% separated</td></tr></table>
      <p class="note">An ordinary collision does not populate the protected circulation; its marked pairs can separate when their marker returns.</p></section>`;
    $("explain").innerHTML = esc(doc.summary) + " The pictures show the probability of constituent A over the plane.";
    $("meta").innerHTML = `<span>carrier</span><b>${doc.size} × ${doc.size} addresses</b><span>ages checked</span><b>${doc.ages.join(", ")}</b><span>contact phase</span><b>${doc.phi.toFixed(6)}</b>` + provenance(doc);
  },
  draw(doc) {
    const dpr = fit(c2d), ctx = c2d.getContext("2d"), W = c2d.width / dpr, H = c2d.height / dpr, N = doc.N;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.fillStyle = "#000"; ctx.fillRect(0, 0, W, H);
    const f = Math.min(doc.F - 1, Math.floor(state.t)), P = Math.min((W - 90) / 2, H * 0.48), y0 = 56;
    const panel = (ch, x0, color, label) => {
      const cell = P / N;
      for (let x = 0; x < N; x++) for (let y = 0; y < N; y++) {
        const v = doc.values[((f * 2 + ch) * N + x) * N + y] / 255;
        ctx.fillStyle = `rgba(${color},${(0.06 + 0.94 * v).toFixed(3)})`; ctx.fillRect(x0 + x * cell + 1, y0 + (N - 1 - y) * cell + 1, cell - 2, cell - 2);
      }
      ctx.strokeStyle = "#1f2c37"; ctx.strokeRect(x0, y0, P, P);
      text(ctx, label, x0, y0 - 10, COL.muted, 12);
      return cell;
    };
    const cell = panel(0, 30, "87,212,154", "the OPH law: the circulation stays whole");
    panel(1, 60 + P, "255,107,107", "independent transport: the pair comes apart");
    const cx = 30 + 4 * cell, cy = y0 + (N - 4) * cell, sgn = f % 2 ? -1 : 1;
    const pts = [[cx, cy], [cx + cell, cy], [cx + cell, cy - cell], [cx, cy - cell]];
    pts.forEach((p, k) => { const q = pts[(k + 1) % 4]; arrow(ctx, p[0] + (q[0] - p[0]) * 0.2, p[1] + (q[1] - p[1]) * 0.2, p[0] + (q[0] - p[0]) * 0.8, p[1] + (q[1] - p[1]) * 0.8, k % 2 ? "#ffb04a" : "#4af1f2", 2.5); });
    dial(ctx, 30 + P - 40, y0 + 40, 22, f * (Math.PI + doc.phi / 2), COL.green, "phase per update π + φ/2");
    const box = [70, y0 + P + 50, W - 110, H - (y0 + P + 50) - 50];
    const { X, Y } = axes(ctx, box, [0, doc.updates], [0, 1.05], "update", "matched weight", { yticks: [0, 0.25, 0.5, 0.75, 1], xticks: ticks(0, doc.updates, 6), font: 12 });
    for (const [arr, col] of [[doc.contact_matched, "#57d49a"], [doc.control_matched, "#ff6b6b"]]) {
      ctx.strokeStyle = col; ctx.lineWidth = 2.5; ctx.beginPath();
      arr.forEach((v, n) => { if (n > f) return; n ? ctx.lineTo(X(n), Y(v)) : ctx.moveTo(X(n), Y(v)); }); ctx.stroke();
    }
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    $("hud").innerHTML = `<b>${esc(doc.step)} · ${esc(doc.title)}</b>`;
    $("clock").innerHTML = `update <b>${f}</b> <small>of ${doc.updates}</small>`;
  },
};

// ---------- OpenWave rows: animated diagrams of each problem ----------
const TAU2 = 2 * Math.PI;
const COL = { ink: "#dce7ee", muted: "#8397a6", line: "#1f2c37", cyan: "#4af1f2", warm: "#ffb04a", violet: "#c38bff", green: "#57d49a", red: "#ff6b6b", blue: "#3a9bff" };
function glow(ctx, x, y, r, color, a = 1) {
  const g = ctx.createRadialGradient(x, y, 0, x, y, r);
  g.addColorStop(0, color); g.addColorStop(1, "rgba(0,0,0,0)");
  ctx.globalAlpha = a; ctx.fillStyle = g; ctx.beginPath(); ctx.arc(x, y, r, 0, TAU2); ctx.fill(); ctx.globalAlpha = 1;
}
function ball(ctx, x, y, r, color, label) {
  glow(ctx, x, y, r * 2.4, color, 0.35);
  ctx.fillStyle = color; ctx.beginPath(); ctx.arc(x, y, r, 0, TAU2); ctx.fill();
  if (label) { ctx.fillStyle = "#000"; ctx.font = `600 ${Math.round(r * 0.9)}px ${MONO()}`; ctx.textAlign = "center"; ctx.textBaseline = "middle"; ctx.fillText(label, x, y + 1); ctx.textAlign = "left"; ctx.textBaseline = "alphabetic"; }
}
function arrow(ctx, x1, y1, x2, y2, color, w = 2) {
  const a = Math.atan2(y2 - y1, x2 - x1), h = 6 + 2 * w;
  ctx.strokeStyle = color; ctx.fillStyle = color; ctx.lineWidth = w;
  ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2 - h * 0.6 * Math.cos(a), y2 - h * 0.6 * Math.sin(a)); ctx.stroke();
  ctx.beginPath(); ctx.moveTo(x2, y2); ctx.lineTo(x2 - h * Math.cos(a - 0.4), y2 - h * Math.sin(a - 0.4)); ctx.lineTo(x2 - h * Math.cos(a + 0.4), y2 - h * Math.sin(a + 0.4)); ctx.closePath(); ctx.fill();
}
function text(ctx, s, x, y, color = COL.muted, size = 13, align = "left") { ctx.fillStyle = color; ctx.font = `${size}px ${MONO()}`; ctx.textAlign = align; ctx.fillText(s, x, y); ctx.textAlign = "left"; }
function wavy(ctx, x1, y1, x2, y2, color, phase) {
  const L = Math.hypot(x2 - x1, y2 - y1), a = Math.atan2(y2 - y1, x2 - x1);
  ctx.strokeStyle = color; ctx.lineWidth = 2; ctx.beginPath();
  for (let s = 0; s <= L; s += 2) { const o = 6 * Math.sin(s / 6 - phase); const x = x1 + s * Math.cos(a) - o * Math.sin(a), y = y1 + s * Math.sin(a) + o * Math.cos(a); s ? ctx.lineTo(x, y) : ctx.moveTo(x, y); }
  ctx.stroke();
}
function dial(ctx, x, y, r, angle, color, label) {
  ctx.strokeStyle = COL.line; ctx.lineWidth = 3; ctx.beginPath(); ctx.arc(x, y, r, 0, TAU2); ctx.stroke();
  ctx.strokeStyle = color; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(x, y); ctx.lineTo(x + r * 0.85 * Math.sin(angle), y - r * 0.85 * Math.cos(angle)); ctx.stroke();
  if (label) text(ctx, label, x, y + r + 18, COL.muted, 12, "center");
}
function plot(ctx, box, f, xr, yr, color, xl, yl) {
  const [x0, y0, w, h] = box; ctx.strokeStyle = COL.line; ctx.lineWidth = 1; ctx.strokeRect(x0, y0, w, h);
  ctx.strokeStyle = color; ctx.lineWidth = 2.5; ctx.beginPath();
  for (let i = 0; i <= 120; i++) { const u = xr[0] + (xr[1] - xr[0]) * i / 120, v = f(u); const X = x0 + w * i / 120, Y = y0 + h - h * (v - yr[0]) / (yr[1] - yr[0]); i ? ctx.lineTo(X, Math.max(y0, Math.min(y0 + h, Y))) : ctx.moveTo(X, Y); }
  ctx.stroke(); if (xl) text(ctx, xl, x0 + w - 4, y0 + h + 16, COL.muted, 12, "right"); if (yl) text(ctx, yl, x0 + 4, y0 - 6);
  return (u, v) => [x0 + w * (u - xr[0]) / (xr[1] - xr[0]), y0 + h - h * (v - yr[0]) / (yr[1] - yr[0])];
}
const DIAGRAMS = {
  charge(ctx, W, H, t) {
    const cx = W * 0.32, cy = H * 0.5, R = Math.min(W, H) * 0.28;
    for (let k = 0; k < 16; k++) { const a = TAU2 * k / 16; arrow(ctx, cx + 22 * Math.cos(a), cy + 22 * Math.sin(a), cx + R * 0.72 * Math.cos(a), cy + R * 0.72 * Math.sin(a), COL.blue, 1.5); }
    ball(ctx, cx, cy, 16, COL.cyan, "e");
    ctx.setLineDash([5, 5]); ctx.strokeStyle = COL.muted; ctx.beginPath(); ctx.arc(cx, cy, R * 0.9, 0, TAU2); ctx.stroke(); ctx.setLineDash([]);
    const a = (t * 0.8) % TAU2; pip(ctx, cx + R * 0.9 * Math.cos(a), cy + R * 0.9 * Math.sin(a), COL.warm);
    text(ctx, `winding around the loop: Q = ${a < TAU2 - 0.05 ? (a / TAU2).toFixed(2) : "1.00"} → 1`, cx, cy + R + 34, COL.ink, 13, "center");
    const sx = W * 0.75; ball(ctx, sx - 40, cy, 11, "rgba(74,241,242,.35)", "½"); ball(ctx, sx + 40, cy, 11, "rgba(74,241,242,.35)", "½");
    ctx.strokeStyle = COL.red; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(sx - 70, cy - 40); ctx.lineTo(sx + 70, cy + 40); ctx.moveTo(sx + 70, cy - 40); ctx.lineTo(sx - 70, cy + 40); ctx.stroke();
    text(ctx, "e → e/2 + e/2 is impossible", sx, cy + 70, COL.red, 13, "center");
  },
  rest(ctx, W, H, t) {
    const cx = W * 0.4, cy = H * 0.5, p = 1 + 0.08 * Math.sin(t * 3);
    glow(ctx, cx, cy, 90 * p, "rgba(74,241,242,.9)", 0.8); ball(ctx, cx, cy, 14, COL.cyan);
    text(ctx, "a localized state the dynamics selects", cx, cy + 120, COL.ink, 13, "center");
    const bx = W * 0.72; ctx.fillStyle = COL.line; ctx.fillRect(bx, cy - 100, 40, 200); ctx.fillStyle = COL.warm; ctx.fillRect(bx, cy + 100 - 140, 40, 140);
    text(ctx, "E = mc²", bx + 20, cy - 110, COL.warm, 15, "center"); text(ctx, "511 keV", bx + 20, cy + 122, COL.muted, 12, "center");
  },
  clock(ctx, W, H, t) {
    const y = H * 0.45, w = 2.2;
    ball(ctx, W * 0.25, y, 14, COL.cyan); dial(ctx, W * 0.25, y - 70, 30, t * w, COL.cyan, "at rest: ω = mc²/ℏ");
    const x = W * 0.45 + ((t * 40) % (W * 0.4)); ball(ctx, x, y + 60, 14, COL.warm); dial(ctx, x, y - 10, 30, t * w * 0.6, COL.warm, "moving: the clock slows");
    arrow(ctx, x + 22, y + 60, x + 60, y + 60, COL.warm, 2);
  },
  persist(ctx, W, H, t) {
    const c = (t % 6) / 6, y = H * 0.5;
    glow(ctx, W * 0.28, y, 50, "rgba(87,212,154,.95)", 0.9); text(ctx, "persists in free evolution", W * 0.28, y + 110, COL.green, 13, "center");
    glow(ctx, W * 0.72, y, 50 + 160 * c, "rgba(255,107,107,.9)", 0.9 * (1 - c) + 0.08); text(ctx, "spreads or collapses (Derrick escape)", W * 0.72, y + 110, COL.red, 13, "center");
  },
  dipole(ctx, W, H, t) {
    const cx = W * 0.45, cy = H * 0.5;
    for (let k = 1; k <= 4; k++) { ctx.strokeStyle = `rgba(58,155,255,${0.6 - k * 0.1})`; ctx.lineWidth = 1.5; for (const s of [-1, 1]) { ctx.beginPath(); ctx.ellipse(cx + s * 40 * k, cy, 40 * k, 30 * k * 1.2, 0, 0, TAU2); ctx.stroke(); } }
    ctx.strokeStyle = COL.warm; ctx.lineWidth = 3; ctx.beginPath(); ctx.ellipse(cx, cy, 40, 14, 0, 0, TAU2); ctx.stroke();
    for (let k = 0; k < 6; k++) { const a = t * 2 + TAU2 * k / 6; pip(ctx, cx + 40 * Math.cos(a), cy + 14 * Math.sin(a), COL.warm); }
    arrow(ctx, cx, cy, cx, cy - 90, COL.cyan, 3); text(ctx, "μ, g ≈ 2", cx + 10, cy - 92, COL.cyan, 14);
  },
  spin(ctx, W, H, t, row) {
    const cx = W * 0.3, cy = H * 0.5, th = (t * 0.9) % (2 * TAU2);
    ctx.save(); ctx.translate(cx, cy); ctx.rotate(th); ctx.fillStyle = "rgba(74,241,242,.25)"; ctx.strokeStyle = COL.cyan; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(0, -70); ctx.lineTo(28, 40); ctx.lineTo(-28, 40); ctx.closePath(); ctx.fill(); ctx.stroke(); ctx.restore();
    text(ctx, `turned by ${(th * 180 / Math.PI).toFixed(0)}°`, cx, cy + 110, COL.ink, 13, "center");
    const s = Math.cos(th / 2), px = W * 0.58, pw = W * 0.34, py = cy - 70;
    const map = plot(ctx, [px, py, pw, 140], u => Math.cos(u / 2), [0, 2 * TAU2], [-1.1, 1.1], COL.violet, "turn 0° → 720°", "");
    text(ctx, "sign of the spin-½ state", px, py - 24, COL.muted, 12);
    const [mx, my] = map(th, s); pip(ctx, mx, my, COL.warm);
    const [x1] = map(TAU2, 0); ctx.strokeStyle = COL.red; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.moveTo(x1, py); ctx.lineTo(x1, py + 140); ctx.stroke(); ctx.setLineDash([]);
    text(ctx, row.number === 7 ? "a full turn: −1 is asked" : "spin ℏ/2", x1 + 4, py - 8 + 0, COL.red, 12);
  },
  annihilation(ctx, W, H, t) {
    const c = (t % 4) / 4, cx = W * 0.5, cy = H * 0.5;
    if (c < 0.5) { const d = (0.5 - c) * 2 * W * 0.35; ball(ctx, cx - d, cy, 13, COL.cyan, "e⁻"); ball(ctx, cx + d, cy, 13, COL.warm, "e⁺"); }
    else { const d = (c - 0.5) * 2 * W * 0.35; glow(ctx, cx, cy, 60 * (1 - (c - 0.5) * 2) + 4, "rgba(255,255,255,.9)"); wavy(ctx, cx, cy, cx, cy - d * 0.8, COL.violet, t * 8); wavy(ctx, cx, cy, cx, cy + d * 0.8, COL.violet, t * 8); text(ctx, "γ", cx + 12, cy - d * 0.8, COL.violet, 16); text(ctx, "γ", cx + 12, cy + d * 0.8, COL.violet, 16); }
    text(ctx, "pair → radiation", cx, H - 24, COL.muted, 13, "center");
  },
  neutral(ctx, W, H, t) {
    const c = (t % 4) / 4, cx = W * 0.3, cy = H * 0.5, d = c * W * 0.55;
    ball(ctx, cx, cy, 15, COL.warm); ball(ctx, cx + d * 0.6, cy - d * 0.35, 9, COL.cyan, "−");
    ctx.setLineDash([3, 6]); ctx.strokeStyle = "rgba(220,231,238,.5)"; ctx.beginPath(); ctx.moveTo(cx, cy); ctx.lineTo(cx + d, cy + d * 0.25); ctx.stroke(); ctx.setLineDash([]);
    glow(ctx, cx + d, cy + d * 0.25, 14, "rgba(220,231,238,.6)"); text(ctx, "ν: no charge, passes through", cx + d, cy + d * 0.25 + 26, COL.ink, 12, "center");
    for (let k = 0; k < 4; k++) { ctx.fillStyle = "rgba(131,151,166,.25)"; ctx.fillRect(W * (0.55 + 0.1 * k), cy - 20, 8, 120); }
  },
  oscillation(ctx, W, H, t) {
    const box = [W * 0.08, H * 0.18, W * 0.84, H * 0.5], L = 4 * Math.PI;
    const f = [u => 1 - 0.8 * Math.sin(u / 2) ** 2, u => 0.55 * Math.sin(u / 2) ** 2 + 0.1 * Math.sin(u) ** 2, u => 0.25 * Math.sin(u / 2) ** 2];
    [COL.cyan, COL.warm, COL.violet].forEach((c, k) => plot(ctx, box, f[k], [0, L], [0, 1.05], c, k ? "" : "distance travelled", k ? "" : "probability of each flavour"));
    const u = (t * 0.8) % L, x = box[0] + box[2] * u / L, p = f.map(g => g(u)), k = p.indexOf(Math.max(...p));
    ball(ctx, x, H * 0.82, 12, [COL.cyan, COL.warm, COL.violet][k], ["e", "μ", "τ"][k]);
  },
  spectrum(ctx, W, H, t, row) {
    const data = row.number === 15 ? [["p", 938.272], ["n", 939.565]] : [["e", 1], ["μ", 206.77], ["τ", 3477.2]];
    const log = row.number !== 15, max = Math.max(...data.map(d => d[1])), x0 = W * 0.15, bw = W * 0.7 / data.length, base = H * 0.82, hmax = H * 0.6, grow = Math.min(1, (t % 6) / 2);
    data.forEach(([n, v], k) => { const frac = log ? (Math.log10(v) + 0.3) / (Math.log10(max) + 0.3) : 0.9 + 0.1 * (v - 938) / 2; const hgt = hmax * frac * grow; ctx.fillStyle = [COL.cyan, COL.warm, COL.violet][k]; ctx.fillRect(x0 + k * bw + bw * 0.2, base - hgt, bw * 0.6, hgt); text(ctx, n, x0 + k * bw + bw * 0.5, base + 22, COL.ink, 16, "center"); text(ctx, log ? `${v}` : `${v} MeV`, x0 + k * bw + bw * 0.5, base - hgt - 8, COL.muted, 12, "center"); });
    text(ctx, log ? "mass ratios, log scale: 1 : 207 : 3477" : "nucleon masses", W * 0.5, H * 0.1, COL.ink, 13, "center");
  },
  quark(ctx, W, H, t) {
    const cx = W * 0.4, cy = H * 0.52, R = Math.min(W, H) * 0.22, cols = ["#ff6b6b", "#57d49a", "#3a9bff"], q = ["+⅔", "+⅔", "−⅓"];
    const pts = [0, 1, 2].map(k => { const a = -Math.PI / 2 + TAU2 * k / 3 + 0.08 * Math.sin(t * 2 + k); return [cx + R * Math.cos(a), cy + R * Math.sin(a)]; });
    pts.forEach((p, k) => { const n = pts[(k + 1) % 3]; ctx.strokeStyle = "rgba(255,176,74,.6)"; ctx.lineWidth = 3 + Math.sin(t * 4 + k); ctx.beginPath(); ctx.moveTo(...p); ctx.lineTo(...n); ctx.stroke(); });
    pts.forEach((p, k) => { ball(ctx, p[0], p[1], 15, cols[k]); text(ctx, q[k], p[0], p[1] - 24, COL.ink, 14, "center"); });
    text(ctx, "charge as a partial winding: ⅔ + ⅔ − ⅓ = 1", W * 0.72, cy, COL.ink, 13, "center");
  },
  baryon(ctx, W, H, t) {
    const cx = W * 0.32, cy = H * 0.5, R = 46, cols = ["#ff6b6b", "#57d49a", "#3a9bff"];
    glow(ctx, cx, cy, 110, "rgba(255,176,74,.25)");
    [0, 1, 2].forEach(k => { const a = t * 1.3 + TAU2 * k / 3; ball(ctx, cx + R * Math.cos(a), cy + R * Math.sin(a), 13, cols[k]); });
    text(ctx, "three colours bound into a colourless state", cx, cy + 130, COL.ink, 12, "center");
    const lx = W * 0.62, lw = W * 0.3, y = e => H * 0.85 - (e - 0.6) * H * 0.7;
    [["three free quarks", 1.0, COL.muted], ["two + one", 0.88, COL.warm], ["bound trimer", 0.7, COL.green]].forEach(([n, e, c]) => { ctx.strokeStyle = c; ctx.lineWidth = 3; ctx.beginPath(); ctx.moveTo(lx, y(e)); ctx.lineTo(lx + lw, y(e)); ctx.stroke(); text(ctx, n, lx, y(e) - 8, c, 12); });
  },
  deuteron(ctx, W, H, t) {
    const cx = W * 0.35, cy = H * 0.5, d = 34 + 4 * Math.sin(t * 2);
    ctx.strokeStyle = "rgba(195,139,255,.5)"; ctx.setLineDash([4, 4]); ctx.beginPath(); ctx.ellipse(cx, cy, 90, 60, 0, 0, TAU2); ctx.stroke(); ctx.setLineDash([]);
    ball(ctx, cx - d, cy, 20, COL.red, "p"); ball(ctx, cx + d, cy, 20, COL.blue, "n");
    text(ctx, "elongated: an electric quadrupole", cx, cy + 92, COL.violet, 12, "center");
    text(ctx, "m_d < m_p + m_n", W * 0.75, cy - 10, COL.ink, 16, "center"); text(ctx, "binding energy 2.22 MeV", W * 0.75, cy + 16, COL.muted, 12, "center");
  },
  nuclear(ctx, W, H, t) {
    const f = A => Math.max(0, 15.8 - 17.8 / Math.cbrt(A) - 0.711 * A * A / (A * Math.cbrt(A)) * (A / 2.1) / A * 0 - 0.711 * (A / 2.2) ** 2 / Math.cbrt(A) / A * 1.0 - 23.7 * 0 ) ;
    const be = A => { const Z = A / (1.98 + 0.0155 * A ** (2 / 3)); return (15.75 * A - 17.8 * A ** (2 / 3) - 0.711 * Z * Z / A ** (1 / 3) - 23.7 * (A - 2 * Z) ** 2 / A) / A; };
    const map = plot(ctx, [W * 0.12, H * 0.15, W * 0.76, H * 0.62], be, [4, 240], [0, 9.5], COL.cyan, "mass number A", "binding energy per nucleon (MeV)");
    const A = 4 + ((t * 25) % 236); const [x, y] = map(A, be(A)); pip(ctx, x, y, COL.warm); text(ctx, "iron peak", ...map(56, 9.2), COL.muted, 12);
  },
  meson(ctx, W, H, t) {
    const cx = W * 0.35, cy = H * 0.5, a = t * 1.5, r = 50;
    const p1 = [cx + r * Math.cos(a), cy + r * Math.sin(a)], p2 = [cx - r * Math.cos(a), cy - r * Math.sin(a)];
    ctx.strokeStyle = "rgba(255,176,74,.7)"; ctx.lineWidth = 4; ctx.beginPath(); ctx.moveTo(...p1); ctx.lineTo(...p2); ctx.stroke();
    ball(ctx, ...p1, 14, COL.red, "q"); ball(ctx, ...p2, 14, COL.cyan, "q̄");
    text(ctx, "π, K: quark and antiquark", cx, cy + 100, COL.ink, 12, "center"); text(ctx, "Λ⁰ → p + π⁻", W * 0.75, cy, COL.ink, 16, "center");
  },
  dark(ctx, W, H, t) {
    const cx = W * 0.45, cy = H * 0.5;
    glow(ctx, cx, cy, 80, "rgba(131,151,166,.8)", 0.6); ball(ctx, cx, cy, 10, "#3a4652");
    const a = t * 0.9; ball(ctx, cx + 120 * Math.cos(a), cy + 70 * Math.sin(a), 7, COL.warm);
    for (let k = 0; k < 3; k++) { const y = cy - 60 + 60 * k, x = ((t * 120 + k * 150) % (W + 100)) - 50; wavy(ctx, x - 60, y, x, y, COL.violet, t * 8); }
    text(ctx, "neutral and massive: light passes, matter orbits", cx, H - 24, COL.ink, 13, "center");
  },
  coulomb(ctx, W, H, t) {
    const y = H * 0.42, x1 = W * 0.22, x2 = W * 0.48;
    for (let k = -3; k <= 3; k++) { ctx.strokeStyle = "rgba(58,155,255,.45)"; ctx.lineWidth = 1.5; ctx.beginPath(); ctx.moveTo(x1, y); ctx.quadraticCurveTo((x1 + x2) / 2, y + k * 45, x2, y); ctx.stroke(); }
    ball(ctx, x1, y, 15, COL.warm, "+"); ball(ctx, x2, y, 15, COL.cyan, "−");
    const map = plot(ctx, [W * 0.62, H * 0.15, W * 0.32, H * 0.55], r => 1 / r, [0.25, 4], [0, 4.2], COL.warm, "distance r", "V ∝ 1/r");
    const r = 0.3 + ((t * 0.6) % 3.6); pip(ctx, ...map(r, 1 / r), COL.ink);
  },
  lorentz(ctx, W, H, t) {
    const v = 0.85 * Math.abs(Math.sin(t * 0.4)), g = 1 / Math.sqrt(1 - v * v), L = W * 0.3, y = H * 0.32;
    ctx.fillStyle = "rgba(87,212,154,.6)"; ctx.fillRect(W * 0.1, y, L, 18); text(ctx, "rod at rest", W * 0.1, y - 8, COL.green, 12);
    ctx.fillStyle = "rgba(255,176,74,.7)"; ctx.fillRect(W * 0.1, y + 60, L / g, 18); text(ctx, `moving at ${v.toFixed(2)} c: length / γ`, W * 0.1, y + 52, COL.warm, 12);
    dial(ctx, W * 0.65, H * 0.4, 36, t * 2, COL.green, "rest clock"); dial(ctx, W * 0.85, H * 0.4, 36, t * 2 / g, COL.warm, `moving clock, γ = ${g.toFixed(2)}`);
    ctx.strokeStyle = "rgba(220,231,238,.35)"; ctx.beginPath(); ctx.moveTo(W * 0.5, H * 0.95); ctx.lineTo(W * 0.38, H * 0.7); ctx.moveTo(W * 0.5, H * 0.95); ctx.lineTo(W * 0.62, H * 0.7); ctx.stroke(); text(ctx, "light cone: nothing faster than c", W * 0.5, H * 0.68, COL.muted, 12, "center");
  },
  confinement(ctx, W, H, t) {
    const sep = 40 + 110 * (0.5 + 0.5 * Math.sin(t * 0.8)), cx = W * 0.3, y = H * 0.42;
    ctx.strokeStyle = "rgba(255,176,74,.85)"; ctx.lineWidth = 7; ctx.beginPath(); ctx.moveTo(cx - sep / 2, y); ctx.lineTo(cx + sep / 2, y); ctx.stroke();
    ball(ctx, cx - sep / 2, y, 14, COL.red, "q"); ball(ctx, cx + sep / 2, y, 14, COL.cyan, "q̄"); text(ctx, "the flux tube stretches", cx, y + 50, COL.ink, 12, "center");
    const map = plot(ctx, [W * 0.6, H * 0.15, W * 0.34, H * 0.6], r => 0.9 * r - 0.3 / r, [0.15, 2.5], [-1.5, 2.4], COL.warm, "separation r", "V = σr − a/r, σ ≈ 1 GeV/fm");
    const r = 0.2 + ((sep - 40) / 110) * 2.2; pip(ctx, ...map(r, 0.9 * r - 0.3 / r), COL.ink);
  },
  running(ctx, W, H, t) {
    const map = plot(ctx, [W * 0.12, H * 0.15, W * 0.76, H * 0.62], q => 1 / (1 + 0.6 * Math.log(q)), [1, 60], [0, 1.05], COL.cyan, "probing scale Q", "coupling α(Q)");
    const q = 1 + ((t * 8) % 59); pip(ctx, ...map(q, 1 / (1 + 0.6 * Math.log(q))), COL.warm);
  },
  decay(ctx, W, H, t, row) {
    const c = (t % 4) / 4, cx = W * 0.4, cy = H * 0.5;
    if (c < 0.35) { ball(ctx, cx, cy, 16, COL.warm, "μ"); glow(ctx, cx, cy, 30 + 20 * c, "rgba(255,176,74,.5)"); }
    else { const d = (c - 0.35) * W * 0.5; ball(ctx, cx + d, cy - d * 0.3, 11, COL.cyan, "e"); pip(ctx, cx - d * 0.6, cy - d * 0.6, COL.ink); pip(ctx, cx - d * 0.5, cy + d * 0.7, COL.ink); text(ctx, "ν", cx - d * 0.6 - 14, cy - d * 0.6, COL.ink, 14); text(ctx, "ν̄", cx - d * 0.5 - 14, cy + d * 0.7, COL.ink, 14); }
    text(ctx, "μ → e + ν + ν̄", W * 0.78, cy, COL.ink, 16, "center");
  },
  beta(ctx, W, H, t) {
    const c = (t % 4) / 4, cx = W * 0.35, cy = H * 0.5;
    if (c < 0.35) ball(ctx, cx, cy, 18, COL.blue, "n");
    else { const d = (c - 0.35) * W * 0.45; ball(ctx, cx, cy, 18, COL.red, "p"); ball(ctx, cx + d, cy + d * 0.3, 10, COL.cyan, "e"); pip(ctx, cx - d * 0.8, cy - d * 0.2, COL.ink); text(ctx, "ν̄", cx - d * 0.8 - 14, cy - d * 0.2, COL.ink, 14); }
    arrow(ctx, W * 0.75, cy + 40, W * 0.75, cy - 40, COL.violet, 3); text(ctx, "spin", W * 0.75 + 8, cy - 40, COL.violet, 12);
    text(ctx, "electrons leave against the spin: parity is broken", W * 0.75, cy + 70, COL.ink, 12, "center");
  },
  newton(ctx, W, H, t) {
    const cx = W * 0.42, cy = H * 0.5, a = t * 0.8, r1 = 50, r2 = 110;
    const p1 = [cx + r1 * Math.cos(a), cy + r1 * Math.sin(a) * 0.6], p2 = [cx - r2 * Math.cos(a), cy - r2 * Math.sin(a) * 0.6];
    ball(ctx, ...p1, 22, COL.warm, "M"); ball(ctx, ...p2, 13, COL.cyan, "m");
    arrow(ctx, p2[0], p2[1], p2[0] + (p1[0] - p2[0]) * 0.3, p2[1] + (p1[1] - p2[1]) * 0.3, COL.ink, 2);
    text(ctx, "F = G M m / r²", W * 0.8, cy, COL.ink, 16, "center");
  },
  metric(ctx, W, H, t) {
    const cx = W * 0.45, cy = H * 0.5;
    ball(ctx, cx, cy, 30, COL.warm, "M");
    ctx.strokeStyle = COL.violet; ctx.lineWidth = 2; ctx.beginPath();
    for (let x = 0; x <= W; x += 4) { const d = x - cx, y = cy - 90 + 50 * Math.exp(-(d * d) / (2 * 120 * 120)) * 1.0; x ? ctx.lineTo(x, y) : ctx.moveTo(x, y); } ctx.stroke();
    const p = (t * 140) % W; pip(ctx, p, cy - 90 + 50 * Math.exp(-((p - cx) ** 2) / (2 * 120 * 120)), COL.violet);
    text(ctx, "light bends", 20, cy - 100, COL.violet, 12);
    dial(ctx, cx, cy + 95, 22, t * 1.4, COL.warm, "near mass: slower"); dial(ctx, W * 0.85, cy + 95, 22, t * 2, COL.green, "far away");
  },
  cosmo(ctx, W, H, t) {
    const s = (t % 8) / 8, a = 0.6 + 0.9 * (s + s * s * 0.8), cx = W * 0.35, cy = H * 0.5;
    for (let i = -3; i <= 3; i++) for (let j = -2; j <= 2; j++) glow(ctx, cx + i * 26 * a, cy + j * 26 * a, 6, "rgba(220,231,238,.9)");
    const map = plot(ctx, [W * 0.66, H * 0.2, W * 0.28, H * 0.55], u => 0.6 + 0.9 * (u + u * u * 0.8), [0, 1], [0.5, 2.3], COL.cyan, "time", "scale factor a(t)");
    pip(ctx, ...map(s, a), COL.warm); text(ctx, "accelerating: Λ", W * 0.8, H * 0.86, COL.warm, 13, "center");
  },
  maxwell(ctx, W, H, t) {
    const y = H * 0.5, x0 = W * 0.08, L = W * 0.84;
    ctx.strokeStyle = COL.line; ctx.beginPath(); ctx.moveTo(x0, y); ctx.lineTo(x0 + L, y); ctx.stroke();
    for (let x = 0; x < L; x += 14) { const ph = x / 40 - t * 3, e = 70 * Math.sin(ph), b = 45 * Math.sin(ph); arrow(ctx, x0 + x, y, x0 + x, y - e, "rgba(58,155,255,.8)", 1.5); arrow(ctx, x0 + x, y, x0 + x + b * 0.5, y + b * 0.35, "rgba(255,176,74,.7)", 1.5); }
    text(ctx, "E", x0, y - 80, COL.blue, 14); text(ctx, "B", x0 + 24, y + 50, COL.warm, 14); text(ctx, "transverse, travelling at c →", x0 + L, y + 100, COL.ink, 13, "right");
  },
  kg(ctx, W, H, t) {
    const map = plot(ctx, [W * 0.1, H * 0.12, W * 0.38, H * 0.66], k => Math.sqrt(1 + k * k), [-3, 3], [0, 3.4], COL.cyan, "wave number k", "ω = √(k² + m²)");
    const k = 2.8 * Math.sin(t * 0.5); pip(ctx, ...map(k, Math.sqrt(1 + k * k)), COL.warm);
    const x0 = W * 0.56, w = W * 0.38, y = H * 0.45;
    ctx.strokeStyle = COL.cyan; ctx.lineWidth = 2; ctx.beginPath();
    for (let i = 0; i <= 200; i++) { const x = i / 200 * 8 - 4, env = Math.exp(-x * x / 2), v = env * Math.cos(2.5 * x - Math.sqrt(1 + 6.25) * t * 2); const X = x0 + w * i / 200, Y = y - 60 * v; i ? ctx.lineTo(X, Y) : ctx.moveTo(X, Y); } ctx.stroke();
    text(ctx, "a matter wave packet", x0 + w / 2, y + 90, COL.ink, 12, "center");
  },
  orbital(ctx, W, H, t) {
    [1, 2, 3].forEach((n, k) => {
      const cx = W * (0.18 + 0.24 * k), cy = H * 0.45, R = 50;
      ctx.strokeStyle = COL.line; ctx.beginPath(); ctx.arc(cx, cy, R, 0, TAU2); ctx.stroke();
      ctx.strokeStyle = [COL.cyan, COL.warm, COL.violet][k]; ctx.lineWidth = 2.5; ctx.beginPath();
      for (let i = 0; i <= 200; i++) { const a = TAU2 * i / 200, r = R + 14 * Math.sin(n * a) * Math.cos(t * 2 / n); i ? ctx.lineTo(cx + r * Math.cos(a), cy + r * Math.sin(a)) : ctx.moveTo(cx + r * Math.cos(a), cy + r * Math.sin(a)); } ctx.stroke();
      ball(ctx, cx, cy, 8, COL.red, ""); text(ctx, `n = ${n}`, cx, cy + R + 30, COL.ink, 13, "center"); text(ctx, `E = −13.6/${n * n} eV`, cx, cy + R + 48, COL.muted, 11, "center");
    });
  },
};
function pip(ctx, x, y, color) { ctx.fillStyle = color; ctx.beginPath(); ctx.arc(x, y, 5, 0, TAU2); ctx.fill(); }

function comparisonChart(ctx, box, cmp) {
  const [x0, y0, w, h] = box;
  ctx.fillStyle = COL.ink; ctx.font = `600 13px ${getComputedStyle(document.body).getPropertyValue("--sans")}`; ctx.fillText(cmp.title, x0, y0 - 10);
  if (cmp.kind === "levels") {
    const lo = Math.min(...cmp.levels.map(l => l[1])) - 8, hi = Math.max(...cmp.levels.map(l => l[2])) + 8, Y = e => y0 + h - h * (e - lo) / (hi - lo);
    cmp.levels.forEach(([name, a, b, color], k) => { const x = x0 + 70 + k * (w - 90) / cmp.levels.length, bw = (w - 90) / cmp.levels.length * 0.7; ctx.fillStyle = color; ctx.globalAlpha = 0.85; ctx.fillRect(x, Y(b), bw, Math.max(3, Y(a) - Y(b))); ctx.globalAlpha = 1; text(ctx, name, x + bw / 2, Y(b) - 8, color, 12, "center"); text(ctx, `${a.toFixed(1)}–${b.toFixed(1)}`, x + bw / 2, Y(a) + 16, COL.muted, 11, "center"); });
    text(ctx, "energy", x0, y0 + 12, COL.muted, 11);
    return;
  }
  const all = cmp.groups.flatMap(g => g[1].map(s => s[1])), signed = cmp.signed || Math.min(...all) < 0;
  const vmax = Math.max(...all.map(Math.abs)) * 1.15, Y = v => signed ? y0 + h / 2 - (h / 2) * v / vmax : y0 + h - h * v / vmax;
  const series = [...new Set(cmp.groups.flatMap(g => g[1].map(s => s[0])))], palette = [COL.cyan, COL.warm, COL.violet, COL.green];
  ctx.strokeStyle = COL.line; ctx.beginPath(); ctx.moveTo(x0, Y(0)); ctx.lineTo(x0 + w, Y(0)); ctx.stroke();
  const gw = w / cmp.groups.length;
  cmp.groups.forEach(([gname, vals], gi) => {
    const bw = gw * 0.7 / Math.max(1, vals.length);
    vals.forEach(([sname, v], si) => { const x = x0 + gi * gw + gw * 0.15 + si * bw; ctx.fillStyle = palette[series.indexOf(sname) % palette.length]; const top = Math.min(Y(v), Y(0)), hh = Math.abs(Y(v) - Y(0)); ctx.fillRect(x + 2, top, bw - 4, Math.max(2, hh)); text(ctx, `${+v.toFixed(4)}${cmp.unit && cmp.unit !== "" ? (cmp.unit.startsWith("/") ? cmp.unit : " " + cmp.unit) : ""}`, x + bw / 2, v >= 0 ? top - 5 : top + hh + 13, COL.ink, 11, "center"); });
    text(ctx, gname, x0 + gi * gw + gw / 2, y0 + h + 18, COL.muted, 12, "center");
  });
  let lx = x0 + Math.min(w * 0.55, ctx.measureText(cmp.title).width + 40);
  series.forEach((s, k) => { ctx.fillStyle = palette[k % palette.length]; ctx.fillRect(lx, y0 - 20, 10, 10); text(ctx, s, lx + 14, y0 - 11, COL.ink, 11); ctx.font = `11px ${MONO()}`; lx += ctx.measureText(s).width + 34; });
}

const rowView = {
  footer: true,
  enter(row) {
    $("legend").innerHTML = "";
    const runs = SITE.runs.filter(r => r.openwave_row === row.number);
    const report = row.report ? `<section><h2>From the current report</h2>${row.report.title ? `<p class="note" style="color:var(--ink);margin-top:0"><b>${esc(row.report.title)}</b></p>` : ""}
      ${(row.report.findings || []).map(s => `<p class="note">${esc(s)}</p>`).join("")}
      ${(row.report.tables || []).map(tb => `${tb.heading ? `<p class="note" style="color:var(--ink)">${esc(tb.heading)}</p>` : ""}<div style="overflow-x:auto"><table><tr>${tb.header.map(h => `<th>${esc(h)}</th>`).join("")}</tr>${tb.rows.map(r => `<tr>${r.map(c => `<td>${esc(c)}</td>`).join("")}</tr>`).join("")}</table></div>`).join("")}</section>` : "";
    $("dyn").innerHTML = `<section><h2>The problem</h2><p class="note" style="color:var(--ink);margin-top:0">${esc(row.explain)}</p>
      <p class="note"><b style="color:var(--ink)">OpenWave's test:</b> ${esc(row.test)}</p>
      <p class="note"><a href="${row.spec_url}" target="_blank" rel="noopener">Problem statement on GitHub</a> · <a href="${row.oph_url}" target="_blank" rel="noopener">OPH column on GitHub</a></p></section>
      <section><h2>OPH's result</h2><p><span class="badge st-${row.status}">${esc(row.status_label)}</span></p>
      <p class="note" style="color:var(--ink)">${row.result ? esc(row.result) : "No passing result for this row is on file."}</p>
      ${runs.map(r => `<a class="go" href="${href(r.path)}">Open the replay: ${esc(r.title)}</a>`).join(" ")}</section>${report}`;
    $("law").innerHTML = `<div class="step">OpenWave · Row ${row.number} · ${esc(row.group)} · ${esc(row.regime)}</div><div class="law-title">${esc(row.name)}</div>
      <div class="kv">${rowCounts()}</div>
      <p class="note">A row passes when a registered finite test of it passed on the OPH law. <a href="${href("openwave/")}">All rows</a></p>`;
    $("explain").innerHTML = "The animation illustrates the physics the row asks for. Where results are on file, the chart below it compares them.";
    $("meta").innerHTML = `<span>criteria</span><b>${esc(OPENWAVE.source.repository)} ${esc(OPENWAVE.source.file)} at ${OPENWAVE.source.commit.slice(0, 8)}</b>`;
  },
  draw(row) {
    const dpr = fit(c2d), ctx = c2d.getContext("2d"), W = c2d.width / dpr, H = c2d.height / dpr;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.fillStyle = "#000"; ctx.fillRect(0, 0, W, H);
    const hasChart = !!row.comparison, dh = hasChart ? H * 0.56 : H - 40;
    ctx.save(); ctx.beginPath(); ctx.rect(0, 30, W, dh); ctx.clip(); ctx.translate(0, 30);
    (DIAGRAMS[row.diagram] || (() => {}))(ctx, W, dh, state.t, row); ctx.restore();
    text(ctx, `Row ${row.number} · ${row.name}`, 16, 22, COL.ink, 13);
    if (hasChart) comparisonChart(ctx, [70, 30 + dh + 40, W - 110, H - dh - 120], row.comparison);
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    $("hud").innerHTML = "";
    $("clock").innerHTML = `<small>animation</small> <b>${state.t.toFixed(1)}</b> s`;
  },
};
const rowOf = n => (OPENWAVE.rows || []).find(r => r.number === n);
function rowLinks(row) {
  if (!row) return "";
  return `<p class="note"><a href="${href(row.path)}">Row ${row.number} page</a> · <a href="${row.spec_url}" target="_blank" rel="noopener">Problem statement on GitHub</a> · <a href="${row.oph_url}" target="_blank" rel="noopener">OPH column</a></p>`;
}
function rowCounts() {
  return [["passed", r => r.status === "pass"], ["open, tested", r => r.status === "open" && r.result], ["no result", r => !r.result]]
    .map(([s, f]) => `<span>${s}</span><b>${OPENWAVE.rows.filter(f).length} rows</b>`).join("");
}
function rowList(current) {
  let group = "", html = "";
  for (const r of OPENWAVE.rows) {
    if (r.group !== group) { group = r.group; html += `<div class="ow-group">${esc(group)}</div>`; }
    const cls = (!r.result ? " open" : "") + (r.number === current ? " current" : "");
    html += `<a class="ow-row${cls}" href="${href(r.path)}"><span class="num">${r.number}</span><span class="nm">${esc(r.name)}<small>${esc(r.test)}</small></span><span class="st st-${r.status}">${esc(r.result ? r.status : "no result")}</span></a>`;
  }
  return html;
}
function showRow(n) {
  const row = rowOf(n); if (!row) return;
  state.run = row; state.view = rowView; row.duration = 24; row.rate = 1; state.t = 0;
  $("empty").style.display = "none"; $("gl").style.display = "none"; $("c2d").style.display = "block"; $("html").style.display = "none";
  document.querySelectorAll(".only-field").forEach(s => s.style.display = "none");
  document.body.classList.remove("static");
  rowView.enter(row);
}
function renderOpenWave() {
  if (!OPENWAVE.rows) { $("owButton").style.display = "none"; return; }
  $("owPanel").innerHTML = `<a class="ow-all" href="${href("openwave/")}">All rows on one page</a>
    <div class="ow-group" style="text-transform:none;letter-spacing:0">Greyed rows have no passing result on the OPH law on file.</div>` + rowList(PAGE.row);
  $("owButton").onclick = e => { e.stopPropagation(); $("owPanel").classList.toggle("on"); };
  document.addEventListener("click", e => { if (!e.target.closest(".ow")) $("owPanel").classList.remove("on"); });
}

// ---------- the contents pages: every experiment, every OpenWave row ----------
function showContents(section) {
  state.run = null;
  $("empty").style.display = "none"; $("gl").style.display = "none"; $("c2d").style.display = "none"; $("html").style.display = "block";
  document.querySelectorAll(".only-field").forEach(s => s.style.display = "none");
  document.body.classList.add("static");
  $("legend").innerHTML = ""; $("hud").innerHTML = "";
  if (section === "openwave" && OPENWAVE.rows) {
    const src = OPENWAVE.source;
    $("html").innerHTML = `<div class="contents"><h3>OpenWave rows</h3>
      <p>The 31 rows of the OpenWave comparison of wave models, from <a href="https://github.com/${esc(src.repository)}/blob/${esc(src.commit)}/${esc(src.file)}" target="_blank" rel="noopener">${esc(src.repository)} ${esc(src.file)}</a>. Each row has its own page with the problem, an animated diagram and OPH's result.</p>
      <div class="ow-list">${rowList()}</div></div>`;
    $("law").innerHTML = `<div class="step">OpenWave</div><div class="law-title">All rows</div><div class="kv">${rowCounts()}</div>
      <p class="note">A row passes when a registered finite test of it passed on the OPH law.</p>`;
    $("explain").innerHTML = "Greyed rows have no passing result on the OPH law on file.";
    $("meta").innerHTML = `<span>criteria</span><b>${esc(src.repository)} ${esc(src.file)} at ${src.commit.slice(0, 8)}</b>`;
  } else {
    const list = EXPERIMENTS.map(r => `<a class="exp" href="${href(r.path)}">${r.step ? `<span class="step">${esc(r.step)}</span>` : ""}<b>${esc(r.title)}</b>${r.summary ? `<small>${esc(r.summary)}</small>` : ""}</a>`).join("");
    $("html").innerHTML = `<div class="contents"><h3>Experiments</h3><p>Every experiment the simulator runs on the OPH law, each on its own page.</p>${list}</div>`;
    $("law").innerHTML = `<div class="step">Experiments</div><div class="law-title">All experiments</div>
      <div class="kv"><span>experiments</span><b>${EXPERIMENTS.length}</b></div>
      <p class="note">The replays of OpenWave rows are on their <a href="${href("openwave/")}">rows' pages</a>.</p>`;
    $("explain").innerHTML = "Each page replays one run file. <code>python -m simulator run all</code> recomputes them.";
    $("meta").innerHTML = "";
  }
}

// ---------- view: a table of clock journeys ----------
const table = {
  footer: false,
  async prepare(doc) { doc.duration = 0; },
  enter(doc) {
    $("legend").innerHTML = "";
    $("html").innerHTML = `<h3>${esc(doc.title)}</h3><p>Each clock's advance on the journey divided by its stationary twin's advance over the same home ticks, read at home at reunion.</p>
      <table style="max-width:760px"><tr><th>route</th><th>reunion</th><th>speed</th><th>forecast</th><th>beat</th><th>optical</th></tr>${doc.journeys.map(j =>
        `<tr><td>(${j.route.join(", ")})</td><td>${j.reunion}</td><td>${j.speed.toFixed(5)}</td><td>${j.forecast.toFixed(5)}</td><td>${j.beat_ratio.toFixed(5)}</td><td>${j.optical_ratio.toFixed(5)}</td></tr>`).join("")}</table>`;
    $("explain").innerHTML = ""; $("meta").innerHTML = provenance(doc);
  },
  draw() {},
};

const VIEWS = { field3d: field, journey, race, scan, formation, dispersion, persistence, "clock-journeys": table };

// ---------- selection and loop ----------
function lawPanel(doc) {
  const l = doc.law, lam = l.equal_phase;
  const rails = l.rail_contact_phases && l.rail_contact_phases.length ? `<span>rail contact phases</span><b>${l.rail_contact_phases.map(v => v.toFixed(6)).join(", ")}</b>` : "";
  $("law").innerHTML = `${doc.step ? `<div class="step">${esc(doc.step)}</div>` : ""}<div class="law-title">${esc(l.title)}</div>
    <span class="badge status-${esc(l.status)}">${esc(l.status)}</span><span class="badge" style="color:var(--muted)">${esc(l.access)}</span>
    <div class="kv"><span>equal-slot phase λ</span><b>${lam[0].toFixed(4)} ${lam[1] >= 0 ? "+" : "−"} ${Math.abs(lam[1]).toFixed(4)}i</b>
    <span>encounter phase χ</span><b>${l.encounter_phase.toFixed(4)} rad</b>${rails}</div>
    <p class="note">${esc(l.note)}</p>`;
}
function select(doc) {
  state.run = doc; state.t = 0; state.view = VIEWS[doc.kind];
  $("runPicker").value = pickerValue(doc);
  $("empty").style.display = "none";
  const isField = doc.kind === "field3d";
  $("gl").style.display = isField ? "block" : "none";
  $("c2d").style.display = !isField && doc.kind !== "clock-journeys" ? "block" : "none";
  $("html").style.display = doc.kind === "clock-journeys" ? "block" : "none";
  document.querySelectorAll(".only-field").forEach(s => s.style.display = isField ? "" : "none");
  $("dyn").innerHTML = "";
  document.body.classList.toggle("static", !state.view.footer);
  lawPanel(doc);
  state.view.enter(doc);
  $("time").value = 0;
}
function frame(now) {
  const doc = state.run, dt = Math.min(0.1, (now - state.last) / 1000);
  state.last = now;
  if (doc) {
    if (state.view.footer && state.playing) {
      state.t += dt * doc.rate;
      if (state.t > doc.duration) state.t = 0;
      $("time").value = Math.round(1000 * state.t / doc.duration);
    }
    if (doc.kind === "field3d" && $("spin").checked) state.yaw += dt * 0.12;
    state.view.draw(doc);
  }
  requestAnimationFrame(frame);
}

// ---------- controls ----------
function segment(id, key) {
  $(id).addEventListener("click", e => {
    const b = e.target.closest("button"); if (!b) return;
    state[key] = b.dataset.v;
    $(id).querySelectorAll("button").forEach(x => x.classList.toggle("on", x === b));
  });
}
segment("frameSeg", "frame"); segment("sourceSeg", "source");
$("play").onclick = () => { state.playing = !state.playing; $("play").textContent = state.playing ? "Pause" : "Play"; };
$("time").oninput = e => { if (state.run && state.run.duration) state.t = state.run.duration * e.target.value / 1000; };
// The picker leads to each experiment's page; a run file opened here stays on this page.
const EXPERIMENTS = SITE.runs.filter(r => !r.openwave_row);
function pickerValue(doc) {
  if (doc.local !== undefined) return "local:" + doc.local;
  const i = SITE.runs.findIndex(r => r.slug === doc.slug && !r.openwave_row);
  return i >= 0 ? "page:" + i : "";
}
function renderPicker() {
  const p = $("runPicker"), listing = PAGE.view === "contents" && PAGE.section === "experiments";
  const listed = PAGE.view === "run" && EXPERIMENTS.some(r => r.slug === PAGE.run);
  if (!listed && !listing) p.add(new Option("Experiments", "", true, true));
  p.add(new Option("All experiments", "contents"));
  SITE.runs.forEach((r, i) => { if (!r.openwave_row) p.add(new Option(r.label, "page:" + i)); });
  if (listing) p.value = "contents";
}
$("runPicker").onchange = e => {
  const v = e.target.value;
  if (v.startsWith("local:")) select(state.runs[+v.slice(6)]);
  else if (v === "contents") location.href = href("experiments/");
  else if (v.startsWith("page:")) location.href = href(SITE.runs[+v.slice(5)].path);
};
$("file").onchange = async e => { let last; for (const f of e.target.files) last = await openRun(new Uint8Array(await f.arrayBuffer()), f.name); if (last) select(last); };
let drag = null;
canvas.addEventListener("pointerdown", e => { drag = [e.clientX, e.clientY]; canvas.setPointerCapture(e.pointerId); });
canvas.addEventListener("pointermove", e => {
  if (!drag) return;
  state.yaw -= (e.clientX - drag[0]) * 0.006; state.pitch = Math.max(-1.4, Math.min(1.4, state.pitch + (e.clientY - drag[1]) * 0.006));
  drag = [e.clientX, e.clientY];
});
canvas.addEventListener("pointerup", () => { drag = null; });
canvas.addEventListener("wheel", e => { e.preventDefault(); state.dist = Math.max(0.9, Math.min(5, state.dist * Math.exp(e.deltaY * 0.001))); }, { passive: false });
const stage = $("stage");
stage.addEventListener("dragover", e => { e.preventDefault(); $("drop").style.display = "flex"; });
stage.addEventListener("dragleave", () => { $("drop").style.display = "none"; });
stage.addEventListener("drop", async e => {
  e.preventDefault(); $("drop").style.display = "none"; let last;
  for (const f of e.dataTransfer.files) last = await openRun(new Uint8Array(await f.arrayBuffer()), f.name);
  if (last) select(last);
});
function glitch() {
  const mark = $("mark");
  if (!mark || matchMedia("(prefers-reduced-motion: reduce)").matches) return;
  const px = a => (Math.random() * 2 - 1) * a + "px";
  const burst = () => {
    for (const [k, a] of [["--rx", 2.5], ["--ry", 1], ["--bx", 2.5], ["--by", 1], ["--sx", 1.5]]) mark.style.setProperty(k, px(a));
    mark.classList.add("glitching");
    setTimeout(() => mark.classList.remove("glitching"), 90 + Math.random() * 160);
    setTimeout(burst, 2500 + Math.random() * 6500);
  };
  setTimeout(burst, 1200);
}
glitch();

(async () => {
  // Optional view in the address: #t=1.2&frame=chart&source=continuum&play=0&spin=0&menu=1
  const view = new URLSearchParams(location.hash.slice(1));
  // Addresses of the earlier one-page app, #run=N and #row=N, lead to their pages.
  if (PAGE.home && (view.has("run") || view.has("row"))) {
    const row = view.has("row") ? rowOf(+view.get("row")) : null, run = view.has("run") ? SITE.runs[+view.get("run")] : null;
    const target = row ? row.path : run ? run.path : null;
    view.delete("run"); view.delete("row");
    if (target) { location.replace(href(target) + (view.toString() ? "#" + view : "")); return; }
  }
  renderPicker();
  renderOpenWave();
  for (const s of EMBEDDED) (await openRun(b64(s))).slug = PAGE.run;
  if (PAGE.view === "run" && state.runs.length) select(state.runs[0]);
  if (PAGE.view === "row") showRow(PAGE.row);
  if (PAGE.view === "contents") showContents(PAGE.section);
  const press = (seg, v) => { const b = $(seg).querySelector(`button[data-v="${v}"]`); if (b) b.click(); };
  if (view.get("frame")) press("frameSeg", view.get("frame"));
  if (view.get("source")) press("sourceSeg", view.get("source"));
  if (view.get("play") === "0") $("play").click();
  if (view.get("spin") === "0") $("spin").checked = false;
  if (view.get("menu")) $("owPanel").classList.add("on");
  if (view.get("t") && state.run && state.run.duration) {
    state.t = Math.min(+view.get("t"), state.run.duration);
    $("time").value = Math.round(1000 * state.t / state.run.duration);
  }
  requestAnimationFrame(frame);
})();
