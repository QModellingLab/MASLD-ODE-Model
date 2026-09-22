Leave-one-out over the silymarin target set  (Reviewer 1 Major 7 / Reviewer 2 Comment 4 /
Reviewer 3 Major 5 follow-up; Supplementary Table S16)

Question
--------
Does the reported stage-dependent response depend on any single silymarin target -- in
particular on CYP2E1, whose supporting evidence is the weakest of the eight (hepatic CYP2E1
is elevated in steatohepatitis, Weltman et al. 1998, but the one direct test found it
unchanged by silibinin, Haddad et al. 2011)?

Method
------
Each of the eight targets is omitted in turn and the remaining seven are inhibited at the
same ki = 0.3 used throughout, under the otherwise unchanged protocol (t = 0-300 h, 3001
points, LSODA, AUC of the active-form trajectory). The full eight-target set is reported on
the first row for reference.

Run
---
    python run_loo.py
Runtime: a few seconds. Writes outputs/leave_one_out_summary.xlsx.

Result
------
The steatosis/steatohepatitis advantage ratio exceeds 1 for all three core outputs in all
eight leave-one-out combinations. Omitting CYP2E1 changes all three ratios by less than 1.5%
(2.4433 -> 2.4074, 2.1048 -> 2.1049, 3.0340 -> 3.0505). The effect is therefore a property
of the target set as a whole, not of any one node.
