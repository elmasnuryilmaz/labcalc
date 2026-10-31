// Static Methods page: formulas, conventions and references.
const METHODS_HTML = `
<section class="card prose">
<h2>Methods</h2>
<p>All calculations run in your browser. The numerical results of every formula below are checked in the
project's automated tests against reference values computed independently with SciPy
(<code>scipy.stats</code>, <code>scipy.optimize.curve_fit</code>); the script that produces them is
<code>scripts/make_reference.py</code> in the repository.</p>

<h3>qPCR relative quantification</h3>
<ol>
<li><b>Technical replicates.</b> Wells with the same biological sample and gene are averaged:
  Ct<sub>s,g</sub> = mean of the included wells. Undetermined or missing Ct values are ignored (and flagged).
  Excluded wells are ignored.</li>
<li><b>Reference genes.</b> With R reference genes, the reference Ct of a sample is the arithmetic mean of their
  Ct values, Ct<sub>ref</sub> = (1/R) Σ Ct<sub>s,r</sub>. Because Ct is on a log<sub>2</sub> scale, this equals using
  the geometric mean of the reference-gene quantities.</li>
<li><b>ΔCt</b> = Ct<sub>target</sub> − Ct<sub>ref</sub> for each biological replicate.</li>
<li><b>ΔΔCt</b> = ΔCt − mean(ΔCt of the control group), for each biological replicate (the control group is the
  calibrator, so its replicates scatter around ΔΔCt = 0).</li>
<li><b>Fold change (Livak)</b> = 2<sup>−ΔΔCt</sup> per replicate (Livak &amp; Schmittgen 2001). This assumes an
  amplification efficiency of 100 % (doubling per cycle) for all genes.</li>
<li><b>Efficiency-corrected ratio (Pfaffl)</b>, used when primer efficiencies are entered. With efficiency
  <i>e</i> in percent, the amplification factor is E = 1 + e/100 (100 % gives E = 2). For a sample s,
  <br>ratio = E<sub>target</sub><sup>ΔCt<sub>target</sub></sup> / ( Π<sub>r</sub> E<sub>r</sub><sup>ΔCt<sub>r</sub></sup> )<sup>1/R</sup>,
  with ΔCt<sub>g</sub> = mean Ct<sub>g</sub>(control group) − Ct<sub>g</sub>(s) (Pfaffl 2001; with several reference genes the
  denominator is their geometric mean). Equivalently, log<sub>2</sub> ratio = −(ΔCq<sub>E</sub>(s) − mean ΔCq<sub>E</sub>(control)),
  where the efficiency-weighted ΔCq<sub>E</sub> = log<sub>2</sub>(E<sub>target</sub>)·Ct<sub>target</sub> − (1/R) Σ log<sub>2</sub>(E<sub>r</sub>)·Ct<sub>r</sub>.
  With all E = 2 this reduces exactly to the Livak result.</li>
<li><b>Group summaries.</b> For each gene and group: n, mean, SD (n − 1 denominator) and SEM = SD/√n of ΔCt, ΔΔCt and
  of the per-replicate fold changes. 2<sup>−(mean ΔΔCt)</sup> (the fold change of the mean, i.e. the geometric mean of the
  fold changes) is also reported. The bar chart shows the mean and SD or SEM of the plotted values (fold change, ratio,
  or their log<sub>2</sub>) with the individual biological replicates as dots.</li>
</ol>

<h3>qPCR statistics</h3>
<p>Tests are performed per gene on the ΔCt values of the biological replicates (on the efficiency-weighted ΔCq<sub>E</sub>
when the Pfaffl method is selected), never on fold changes, because ΔCt is on a log scale and closer to normally
distributed. Technical replicates are averaged first and are never treated as independent observations.</p>
<ul>
<li><b>Two groups:</b> Welch's unequal-variance t-test, t = (m<sub>1</sub> − m<sub>2</sub>) / √(s<sub>1</sub>²/n<sub>1</sub> + s<sub>2</sub>²/n<sub>2</sub>),
  with Welch–Satterthwaite degrees of freedom; two-sided p-value.</li>
<li><b>More than two groups:</b> one-way ANOVA (F test), followed by either
  <ul>
  <li><b>Tukey HSD</b> for all pairs (Tukey–Kramer form for unequal group sizes):
    q = |m<sub>i</sub> − m<sub>j</sub>| / √(MSE/2 · (1/n<sub>i</sub> + 1/n<sub>j</sub>)), p from the studentized range distribution with
    k groups and N − k degrees of freedom; or</li>
  <li><b>Dunnett's test</b> for each group against the control: t<sub>i</sub> = (m<sub>i</sub> − m<sub>0</sub>) / √(MSE (1/n<sub>i</sub> + 1/n<sub>0</sub>)),
    two-sided single-step p-value from the joint (multivariate t) distribution of the k − 1 statistics, whose correlation
    λ<sub>i</sub>λ<sub>j</sub> with λ<sub>i</sub> = √(n<sub>i</sub> / (n<sub>i</sub> + n<sub>0</sub>)) reduces the probability to a two-dimensional
    integral evaluated numerically.</li>
  </ul></li>
</ul>
<p>The post-hoc tests use the pooled ANOVA mean square error (MSE). The ANOVA-based tests assume similar variances
across groups. p-values are not adjusted across genes.</p>

<h3>qPCR quality control</h3>
<ul>
<li>SD of technical replicates (Ct) above the "note" threshold (default 0.3 cycles) gives a warning, above the "flag"
  threshold (default 0.5) an error with an offer to exclude the replicate furthest from the median.</li>
<li>Ct above 35 (adjustable) is flagged as late amplification; undetermined or missing Ct values are flagged and ignored.</li>
<li>Any amplification in a no-template control (NTC) is flagged (an error when Ct is below the late-amplification
  threshold). Check the melt curve for primer dimers.</li>
<li>Samples missing a valid reference or target value are dropped from that calculation and listed.</li>
</ul>

<h3>Cell viability (MTT, WST-1)</h3>
<ol>
<li>Optional reference wavelength: A = A<sub>measure</sub> − A<sub>reference</sub>, well by well.</li>
<li>Blank: the mean of the blank (medium only) wells is subtracted from every well.</li>
<li>% viability = 100 · (A<sub>well</sub> − mean A<sub>blank</sub>) / (mean A<sub>untreated control</sub> − mean A<sub>blank</sub>).</li>
<li>Replicate wells at each concentration are summarised as mean and SD; the curve is fitted to the individual wells.</li>
</ol>

<h3>Dose-response fit</h3>
<p>Four-parameter logistic model on x = log<sub>10</sub>(concentration):</p>
<p class="formula">y = Bottom + (Top − Bottom) / (1 + 10<sup>Hill · (x − log<sub>10</sub> IC50)</sup>)</p>
<p>Parameters are estimated by least squares with the Levenberg–Marquardt algorithm, started from several initial values
(the best fit is kept). Results are reported with Hill &gt; 0, so Top is the plateau at low concentration and Bottom the
plateau at high concentration. This is the same model as LL.4 in the R package <i>drc</i>. Standard errors come from the
asymptotic covariance matrix s²(JᵀJ)<sup>−1</sup>, with s² = SSR/(n − p) and p the number of fitted parameters. The 95%
confidence interval of the IC50 is log<sub>10</sub>IC50 ± t<sub>0.975, n−p</sub> · SE, back-transformed, so it is asymmetric
on the concentration scale. R² = 1 − SSR/SST. Top and/or Bottom can be fixed at 100 % and 0 % when the data do not define
a plateau.</p>
<p>The reported IC50 is the <i>relative</i> IC50 (the midpoint between the fitted plateaus). The <i>absolute</i> IC50, the
concentration at which the fitted curve crosses 50 % viability, is reported alongside when the curve crosses 50 %.
Fits are flagged when the IC50 lies outside the tested concentration range (an extrapolation), when the 95% CI spans more
than two orders of magnitude, when the Hill slope is very steep, when the curve does not reach 50 %, or when there are no
residual degrees of freedom.</p>

<h3>References</h3>
<ol class="refs">
<li>Livak KJ, Schmittgen TD (2001). Analysis of relative gene expression data using real-time quantitative PCR and the
  2<sup>−ΔΔCT</sup> method. <i>Methods</i> 25(4):402–408. doi:10.1006/meth.2001.1262</li>
<li>Pfaffl MW (2001). A new mathematical model for relative quantification in real-time RT-PCR.
  <i>Nucleic Acids Research</i> 29(9):e45. doi:10.1093/nar/29.9.e45</li>
<li>Vandesompele J, De Preter K, Pattyn F, et al. (2002). Accurate normalization of real-time quantitative RT-PCR data by
  geometric averaging of multiple internal control genes. <i>Genome Biology</i> 3(7):research0034.
  doi:10.1186/gb-2002-3-7-research0034</li>
<li>Bustin SA, Benes V, Garson JA, et al. (2009). The MIQE guidelines: minimum information for publication of quantitative
  real-time PCR experiments. <i>Clinical Chemistry</i> 55(4):611–622. doi:10.1373/clinchem.2008.112797</li>
<li>Welch BL (1947). The generalization of ‘Student's’ problem when several different population variances are involved.
  <i>Biometrika</i> 34(1–2):28–35. doi:10.1093/biomet/34.1-2.28</li>
<li>Kramer CY (1956). Extension of multiple range tests to group means with unequal numbers of replications.
  <i>Biometrics</i> 12(3):307–310. doi:10.2307/3001469</li>
<li>Dunnett CW (1955). A multiple comparison procedure for comparing several treatments with a control.
  <i>Journal of the American Statistical Association</i> 50(272):1096–1121. doi:10.1080/01621459.1955.10501294</li>
<li>Levenberg K (1944). A method for the solution of certain non-linear problems in least squares.
  <i>Quarterly of Applied Mathematics</i> 2(2):164–168. doi:10.1090/qam/10666</li>
<li>Marquardt DW (1963). An algorithm for least-squares estimation of nonlinear parameters.
  <i>Journal of the Society for Industrial and Applied Mathematics</i> 11(2):431–441. doi:10.1137/0111030</li>
</ol>
<p class="hint">The numerical output of every procedure is validated against SciPy.</p>
</section>`;

const ABOUT_HTML = `
<section class="card prose">
<h2>Privacy</h2>
<p><b>Your data never leave your computer.</b> labcalc is a static web page: files you open are read by your browser
and every calculation runs locally in JavaScript. There is no server, no account, no analytics, no cookies and no
external fonts or scripts. The published page carries a Content-Security-Policy with <code>connect-src 'none'</code>,
which makes the browser itself block any network request the page might try to make. Downloads (CSV, SVG, PNG) are
generated in the browser. Closing the tab discards everything.</p>
<h2>About</h2>
<p>labcalc is free to use. It can also be run offline: build it from the repository and serve the
<code>dist</code> folder with any local static file server (see the README). Example datasets in the app and repository are
<b>synthetic</b> (simulated) and are only meant for trying the tool.</p>
<p>Results are provided without warranty. Check the QC flags and the fit diagnostics before using numbers in a
publication, and report the method used (Livak or Pfaffl, reference genes, efficiencies, statistical test).</p>
</section>`;

export function mountMethods(el: HTMLElement): void {
  el.innerHTML = METHODS_HTML;
}
export function mountAbout(el: HTMLElement): void {
  el.innerHTML = ABOUT_HTML;
}
