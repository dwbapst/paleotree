07-06-24:
Double-checked to ensure there are package anchors for all Rd \link{} targets. None noted in R CHECK locally, on winbuilder, or using Github Actions.

07-05-24: 
One new feature under the hood, passes R CHECK both locally and on win builder.

08-21-22:
- Fixed broken HTML code that was causing an issue in one Rd help file and reported to maintainer on 08-19-22. Checks out without issue on win-builder using R-devel now.

04-21-22:
- Broken URL fixed. Winbuilder also was back up -- finished checks without issue.
- All functions that call on external Internet resources now fail gracefully if the resource is not available when examples are run, including examples in 'donttest' sections. This should address concerns raised in emails from CRAN about paleotree in April 2022. R CHECK passes locally on Ubuntu, and on a remote Windows environment (Winbuilder appears to be down at the moment).

06-04-19: 
Lots of revision to pre-existing code, passes R CHECK locally, and on win builder.

10-01-18: 
Some new features, passes R CHECK both locally and on win builder.
