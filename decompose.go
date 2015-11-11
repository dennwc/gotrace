package gotrace

import (
	"fmt"
	"sync"
)

var detrand_t = [256]byte{ // non-linear sequence: constant term of inverse in GF(8), mod x^8+x^4+x^3+x+1
	0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1,
	0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0,
	0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
	1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1,
	0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0,
	0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0,
	0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0,
	0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1,
	1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0,
	0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1,
	1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
}

// deterministically and efficiently hash (x,y) into a pseudo-random bit
func detrand(x, y int) bool {
	t := detrand_t[:]
	z := ((0x04b3e375 * x) ^ y) * 0x05a8ef93
	z = int(t[z&0xff]) ^ int(t[(z>>8)&0xff]) ^ int(t[(z>>16)&0xff]) ^ int(t[(z>>24)&0xff])
	return z != 0
}

// set the excess padding to 0
func (bm *Bitmap) clearexcess() {
	if Word(bm.W)%wordBits != 0 {
		mask := allBits << (wordBits - (Word(bm.W) % wordBits))
		for y := 0; y < bm.H; y++ {
			*(bm.index(bm.W, y)) &= mask
		}
	}
}

// return the "majority" value of bitmap bm at intersection (x,y).
// We assume that the bitmap is balanced at "radius" 1.
func (bm *Bitmap) majority(x, y int) bool {
	for i := 2; i < 5; i++ { // check at "radius" i
		ct := 0
		for a := -i + 1; a <= i-1; a++ {
			if bm.Get(x+a, y+i-1) {
				ct++
			} else {
				ct--
			}
			if bm.Get(x+i-1, y+a-1) {
				ct++
			} else {
				ct--
			}
			if bm.Get(x+a-1, y-i) {
				ct++
			} else {
				ct--
			}
			if bm.Get(x-i, y+a) {
				ct++
			} else {
				ct--
			}
		}
		if ct > 0 {
			return true
		} else if ct < 0 {
			return false
		}
	}
	return false
}

// decompose image into paths

// efficiently invert bits [x,infty) and [xa,infty) in line y.
// Here xa must be a multiple of wordBits.
func (bm *Bitmap) xor_to_ref(x, y, xa int) {
	xhi := x & -int(wordBits)
	xlo := x & int(wordBits-1) // = x % BM_WORDBITS

	if xhi < xa {
		for i := xhi; i < xa; i += int(wordBits) {
			*(bm.index(i, y)) ^= allBits
		}
	} else {
		for i := xa; i < xhi; i += int(wordBits) {
			*(bm.index(i, y)) ^= allBits
		}
	}
	// note: the following "if" is needed because x86 treats a<<b as
	// a<<(b&31). I spent hours looking for this bug.
	if xlo != 0 {
		*bm.index(xhi, y) ^= (allBits << (wordBits - Word(xlo)))
	}
}

//	a path is represented as an array of points, which are thought to
//	lie on the corners of pixels (not on their centers). The path point
//	(x,y) is the lower left corner of the pixel (x,y). Paths are
//	represented by the len/pt components of a path_t object (which
//	also stores other information about the path)

// xor the given pixmap with the interior of the given path.
// Note: the path must be within the dimensions of the pixmap.
func (bm *Bitmap) xor_path(p *Path) {
	if len(p.priv.Pt) <= 0 { // a path of length 0 is silly, but legal
		return
	}

	y1 := int(p.priv.Pt[len(p.priv.Pt)-1].Y)

	xa := int(p.priv.Pt[0].X) & -int(wordBits)
	for k := 0; k < len(p.priv.Pt); k++ {
		x := int(p.priv.Pt[k].X)
		y := int(p.priv.Pt[k].Y)

		if y != y1 {
			// efficiently invert the rectangle [x,xa] x [y,y1]
			bm.xor_to_ref(x, min(y, y1), xa)
			y1 = y
		}
	}
}

var ptPool = &sync.Pool{}

//	compute a path in the given pixmap, separating black from white.
//	Start path at the point (x0,x1), which must be an upper left corner
//	of the path. Also compute the area enclosed by the path. Return a
//	new path_t object, or NULL on error (note that a legitimate path
//	cannot have length 0). Sign is required for correct interpretation
//	of turnpolicies.
func (bm *Bitmap) findpath(x0, y0 int, sign int, turnpolicy TurnPolicy) (*Path, error) {
	var (
		xn0, yn0   = int32(x0), int32(y0)
		x, y       = xn0, yn0
		dirx, diry = int32(0), int32(-1)
		area       int
		pt         [][2]int32
	)
	if arr := ptPool.Get(); arr != nil {
		pt = arr.([][2]int32)
		pt = pt[0:0]
	}
	defer func() {
		if len(pt) > 1<<16 {
			ptPool.Put(pt)
		}
	}()
	const limit = 1 << 24
	for i := 0; ; i++ {
		if i >= limit {
			return nil, fmt.Errorf("limit reached")
		}
		// add point to path
		pt = append(pt, [2]int32{int32(x), int32(y)})

		// move to next point
		x += dirx
		y += diry
		area += int(x) * int(diry)

		// path complete?
		if x == xn0 && y == yn0 {
			break
		}

		// determine next direction
		c := bm.Get(int(x+(dirx+diry-1)/2), int(y+(diry-dirx-1)/2))
		d := bm.Get(int(x+(dirx-diry-1)/2), int(y+(diry+dirx-1)/2))

		if c && !d { // ambiguous turn
			if turnpolicy == TurnRight ||
				(turnpolicy == TurnBlack && sign == +1) ||
				(turnpolicy == TurnWhite && sign == -1) ||
				(turnpolicy == TurnRandom && detrand(int(x), int(y))) ||
				(turnpolicy == TurnMajority && bm.majority(int(x), int(y))) ||
				(turnpolicy == TurnMinority && !bm.majority(int(x), int(y))) {
				dirx, diry = diry, -dirx // right turn
			} else {
				dirx, diry = -diry, dirx // left turn
			}
		} else if c { // right turn
			dirx, diry = diry, -dirx
		} else if !d { // left turn
			dirx, diry = -diry, dirx
		}
	} // while this path
	path := make([]point, len(pt))
	for i, p := range pt {
		path[i] = point{int(p[0]), int(p[1])}
	}
	return &Path{
		priv: &privPath{Pt: path},
		Area: area, Sign: sign,
	}, nil
}

//	find the next set pixel in a row <= y. Pixels are searched first
//	left-to-right, then top-down. In other words, (x,y)<(x',y') if y>y'
//	or y=y' and x<x'. If found, return 0 and store pixel in
//	(*xp,*yp). Else return 1. Note that this function assumes that
//	excess bytes have been cleared with bm_clearexcess.
func (bm *Bitmap) findnext(xp, yp *int) bool {
	x0 := (*xp) & ^int(wordBits-1)
	for y := *yp; y >= 0; y-- {
		for x := x0; x < bm.W; x += int(wordBits) {
			if *bm.index(x, y) != 0 {
				for !bm.Get(x, y) {
					x++
				}
				// found
				*xp = x
				*yp = y
				return true
			}
		}
		x0 = 0
	}
	// not found
	return false
}

func (p *Path) clearPriv() {
	p.priv = nil
}

//	Decompose the given bitmap into paths. Returns a linked list of
//	path_t objects with the fields len, pt, area, sign filled
//	in. Returns 0 on success with plistp set, or -1 on error with errno
//	set.
func (bm *Bitmap) toPathList(param *Params) (out []Path, err error) {
	bm1 := bm.Clone()
	// be sure the byte padding on the right is set to 0, as the fast pixel search below relies on it
	bm1.clearexcess()
	// iterate through components
	x, y := 0, bm1.H-1

	type par struct {
		x, y int
	}

	var sign int
	mp := make(map[par]int)
	for bm1.findnext(&x, &y) {
		// calculate the sign by looking at the original
		if bm.Get(x, y) {
			sign = +1
		} else {
			sign = -1
		}
		// calculate the path
		if c := mp[par{x, y}]; c > 5 {
			return out, fmt.Errorf("killed on endless loop")
		}
		p, err := bm1.findpath(x, y+1, sign, param.TurnPolicy)
		if err != nil {
			return out, err
		}
		cc := mp[par{x, y}]
		cc++
		mp[par{x, y}] = cc
		// update buffered image
		bm1.xor_path(p)
		// if it's a turd, eliminate it, else append it to the list
		if p.Area > param.TurdSize {
			processOnePath(p, param)
			out = append(out, *p)
		}
	}

	out = bm1.pathlist_to_tree(out)

	var clear func(p *Path)
	clear = func(p *Path) {
		p.clearPriv()
		for i := range p.Childs {
			clear(&p.Childs[i])
		}
	}
	for i := range out {
		clear(&out[i])
	}
	return
}

type bbox struct {
	x0, x1, y0, y1 int
}

// Find the bounding box of a given path
func setbbox_path(p *Path) (bbox bbox) {
	//	if len(p.priv.Pt) == 0 {
	//		return
	//	}
	bbox.y0 = p.priv.Pt[0].Y
	bbox.y1 = 0
	bbox.x0 = p.priv.Pt[0].X
	bbox.x1 = 0

	for k := range p.priv.Pt {
		x := p.priv.Pt[k].X
		y := p.priv.Pt[k].Y
		if x < bbox.x0 {
			bbox.x0 = x
		}
		if x > bbox.x1 {
			bbox.x1 = x
		}
		if y < bbox.y0 {
			bbox.y0 = y
		}
		if y > bbox.y1 {
			bbox.y1 = y
		}
	}
	return
}

// clear the bm, assuming the bounding box is set correctly (faster than clearing the whole bitmap)
func (bm *Bitmap) clear_bm_with_bbox(bbox bbox) {
	imin := (bbox.x0 / int(wordBits))
	imax := ((bbox.x1 + int(wordBits) - 1) / int(wordBits))

	for y := bbox.y0; y < bbox.y1; y++ {
		for i := imin; i < imax; i++ {
			bm.scanline(y)[i] = 0
		}
	}
}

//	Give a tree structure to the given path list, based on "insideness"
//	testing. I.e., path A is considered "below" path B if it is inside
//	path B. The input pathlist is assumed to be ordered so that "outer"
//	paths occur before "inner" paths. The tree structure is stored in
//	the "childlist" and "sibling" components of the path_t
//	structure. The linked list structure is also changed so that
//	negative path components are listed immediately after their
//	positive parent.  Note: some backends may ignore the tree
//	structure, others may use it e.g. to group path components. We
//	assume that in the input, point 0 of each path is an "upper left"
//	corner of the path, as returned by bm_to_pathlist. This makes it
//	easy to find an "interior" point. The bm argument should be a
//	bitmap of the correct size (large enough to hold all the paths),
//	and will be used as scratch space. Return 0 on success or -1 on
//	error with errno set.

func (bm *Bitmap) pathlist_to_tree(inp []Path) []Path {
	type path struct {
		Path      *Path
		next      *path
		childlist *path
		sibling   *path
	}
	var (
		cur, head, p *path
	)
	bm.Clear(false)

	var plist *path
	{
		cur := &plist
		for i := range inp {
			lp := &path{Path: &inp[i]}
			*cur = lp
			cur = &lp.next
		}
	}

	// save original "next" pointers
	for cur = plist; cur != nil; cur = cur.next {
		cur.sibling = cur.next
	}

	heap := plist

	//	the heap holds a list of lists of paths. Use "childlist" field
	//	for outer list, "next" field for inner list. Each of the sublists
	//	is to be turned into a tree. This code is messy, but it is
	//	actually fast. Each path is rendered exactly once. We use the
	//	heap to get a tail recursive algorithm: the heap holds a list of
	//	pathlists which still need to be transformed.

	list_insert_beforehook := func(elt **path, hook ***path) {
		(*elt).next = **hook
		**hook = *elt
		*hook = &((*elt).next)
	}
	list_append := func(list, elt **path) {
		var hook **path
		for hook = list; *hook != nil; hook = &((*hook).next) {
		}
		(*elt).next = *hook
		*hook = *elt
	}
	list_forall_unlink := func(elt, list **path, fnc func(elt *path) bool) {
		chk := func(elt *path) bool {
			if elt != nil {
				*list, elt.next = elt.next, nil
				return true
			}
			return false
		}
		for *elt = *list; chk(*elt); *elt = *list {
			if !fnc(*elt) {
				break
			}
		}
	}

	for heap != nil {
		// unlink first sublist
		cur = heap
		heap = heap.childlist
		cur.childlist = nil

		// unlink first path
		head = cur
		cur = cur.next
		head.next = nil

		// render path
		bm.xor_path(head.Path)
		bbox := setbbox_path(head.Path)

		// now do insideness test for each element of cur; append it to
		// head.childlist if it's inside head, else append it to
		// head.next.
		hook_in := &(head.childlist)
		hook_out := &(head.next)

		// list_forall_unlink(p, cur)
		list_forall_unlink(&p, &cur, func(p *path) bool {
			if p.Path.priv.Pt[0].Y <= bbox.y0 {
				list_insert_beforehook(&p, &hook_out)
				// append the remainder of the list to hook_out
				*hook_out = cur
				return false
			}
			if bm.Get(p.Path.priv.Pt[0].X, p.Path.priv.Pt[0].Y-1) {
				list_insert_beforehook(&p, &hook_in)
			} else {
				list_insert_beforehook(&p, &hook_out)
			}
			return true
		})

		// clear bm
		bm.clear_bm_with_bbox(bbox)

		// now schedule head.childlist and head.next for further processing
		if head.next != nil {
			head.next.childlist = heap
			heap = head.next
		}
		if head.childlist != nil {
			head.childlist.childlist = heap
			heap = head.childlist
		}
	}

	// copy sibling structure from "next" to "sibling" component
	for p = plist; p != nil; p, p.sibling = p.sibling, p.next {
	}

	//	reconstruct a new linked list ("next") structure from tree
	//	("childlist", "sibling") structure. This code is slightly messy,
	//	because we use a heap to make it tail recursive: the heap
	//	contains a list of childlists which still need to be
	//	processed.
	heap = plist
	if heap != nil {
		heap.next = nil // heap is a linked list of childlists
	}
	plist = nil
	plist_hook := &plist
	for heap != nil {
		heap1 := heap.next
		for p = heap; p != nil; p = p.sibling {
			// p is a positive path
			// append to linked list
			list_insert_beforehook(&p, &plist_hook)

			// go through its children
			for p1 := p.childlist; p1 != nil; p1 = p1.sibling {
				// append to linked list
				list_insert_beforehook(&p1, &plist_hook)
				// append its childlist to heap, if non-empty
				if p1.childlist != nil {
					//list_append(path_t, heap1, p1.childlist)
					list_append(&heap1, &p1.childlist)
				}
			}
		}
		heap = heap1
	}

	// Now, form Go slice structure from that mess
	var convPaths func(cur *path) []Path
	convPaths = func(cur *path) (out []Path) {
		if cur == nil {
			return nil
		}
		for ; cur != nil; cur = cur.sibling {
			cur.Path.Childs = convPaths(cur.childlist)
			out = append(out, *cur.Path)
		}
		return
	}
	return convPaths(plist)
}

func (bm *Bitmap) print() {
	var w Word
	for i := range bm.Map {
		w += bm.Map[i]
	}
	fmt.Printf("bm: %d\n", w)
}
