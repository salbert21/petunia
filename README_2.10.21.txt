NOTE: figured out today that Lu uses a heading angle of due-North for 0, whereas I have been using a heading angle of
	due-East for 0 (both measured clockwise). Unfortunately the code is currently a mix of both conventions.
	Code from Lu or ASEN 6015, including: cart2sph, sph2cart, vr2viCart, and sphericalEntryEOMs all use the 0-North
	convention.

TODO: In the long-term, code badly needs cleanup and to have a single uniform convention. In addition to heading-angle
	difference, Lu's EOMs also require radians and m, m/s where I have been using km, km/s, deg.

      In the short-term, I need to finish this FNPAG/FNPAG work before doing the above. I will leave both sets of code
	active, even though they are inconsistent. For now, simply adding hda + 90deg to the spherical state when going
	from my state to a state for Lu's EOMs (or vice-versa) solves the problem. Check out test_conversions.py and
	test_sphericalEOMs.py for basic comparison checks between algorithms. Delete these scripts only once code is
	uniform and consistent.

NOTE 2: also realized that some code was modifying the value of attributes of a Param instance inside a function. Like
	lists, this user-defind class object IS MUTABLE, and thus is BEING PASSED BY REFERENCE! Therefore, short of a
	significant refactor, this requires a new rule for Petunia: VALUES IN A PARAM INSTANCE SHOULD NOT BE CHANGED
	DURING A SIMULATION. THIS IS ONLY FOR CONTAINING CONSTANTS (vehicle params, tolerances, initial state, etc.),
	NOT for containing any current state information.