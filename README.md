# Causal Meta-Mediation Analysis
## KDD 2020: Accepted and Oral Presentation (Selection Rate: 5.8%, 44 out of 756)

## Citation
Zenan Wang, Xuan Yin, Tianbo Li, and Liangjie Hong. 2020. Causal Meta-Mediation Analysis: Inferring Dose-Response Function From Summary Statistics of Many Randomized Experiments. In *Proceedings of the 26th ACM SIGKDD Conference on Knowledge Discovery and Data Mining (KDD ’20), August 23-27, 2020, Virtual Event, CA, USA*. ACM, New York, NY, USA, 11 pages. https://doi.org/10.1145/3394486.3403313

## <ins>[The Author Version of the Paper](https://github.com/znwang25/cmma/blob/master/KDD2020_CMMA_author_version.pdf)</ins>

## Implementation
[link]()

## Abstract
It is common in the internet industry to use offline-developed algorithms to power online products that contribute to the success of a business.  Offline-developed algorithms are guided by offline evaluation metrics, which are often different from online business key performance indicators (KPIs).  To maximize business KPIs, it is important to pick a north star among all available offline evaluation metrics.  By noting that online products can be measured by online evaluation metrics, the online counterparts of offline evaluation metrics, we decompose the problem into two parts.  As the offline A/B test literature works out the first part: counterfactual estimators of offline evaluation metrics that move the same way as their online counterparts, we focus on the second part: causal effects of online evaluation metrics on business KPIs.  The north star of offline evaluation metrics should be the one whose online counterpart causes the most significant lift in the business KPI.  We model the online evaluation metric as a mediator and formalize its causality with the business KPI as dose-response function (DRF).  Our novel approach, causal meta-mediation analysis, leverages summary statistics of many existing randomized experiments to identify, estimate, and test the mediator DRF.  It is easy to implement and to scale up, and has many advantages over the literature of mediation analysis and meta-analysis.  We demonstrate its effectiveness by simulation and implementation on real data.

## Keywords
causal inference; meta-analysis; mediation analysis; experiment; dose-response function; A/B test; evaluation metric; business KPI
