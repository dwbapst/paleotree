#####

04-21-22:
-Broken URL fixed. Winbuilder also was back up -- finished checks without issue.
-All functions that call on external Internet resources now fail gracefully if the resource is not available when examples are run, including examples in 'donttest' sections. This should address concerns raised in emails from CRAN about paleotree in April 2022. R CHECK passes locally on Ubuntu, and on a remote Windows environment (Winbuilder appears to be down at the moment).

06-04-19: 
Lots of revision to pre-existing code, passes R CHECK locally, and on win builder.

10-01-18: 
Some new features, passes R CHECK both locally and on win builder.
