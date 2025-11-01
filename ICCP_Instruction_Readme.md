### **What These Files Do (A Simple Guide)**

We’re trying to take a **"messed-up route"** (a path from a "dumb" navigation system, the INS) and figure out where it actually belongs on a **"real map"** (the geomagnetic map). We use the magnetic readings from the route to help us find the right spot.

### 1. **Main_ICCP_Geomagnetic.m**

* **(Job: The "Boss")**

**Purpose**: This is the only file you need to run. It's the "boss" that tells all the other files what to do and in what order.

**What it does**:

* Calls the **"Map Maker"** (ICCP_Setup_Simulation) to get a map and a messed-up route to fix.
* Hands off that map and route to the **"Problem Solver"** (ICCP_Algorithm) to fix it.
* Creates the final graphs showing:

  * The **"messed-up"** route
  * The **"real"** route
  * The **"fixed"** route

**Relationships**:

* **Calls**: ICCP_Setup_Simulation.m, ICCP_Algorithm.m
* **Called By**: (You, the user)

---

### 2. **ICCP_Setup_Simulation.m**

* **(Job: The "Map Maker" / Puzzle Creator)**

**Purpose**: This file’s job is to create the puzzle the algorithm will solve.

**What it does**:

* Builds the **"Real Map"**: Generates a geomagnetic map with peaks to simulate hills and valleys.
* Creates the **"Real Route"**: Draws the perfect path (the "ground truth").
* Creates the **"Messed-up Route"**: Distorts the perfect path by rotating, moving, and stretching it.
* Collects **"Clues"**: Records the magnetic readings at each point on the messed-up route.

**Relationships**:

* **Calls**: MATLAB built-in functions (like `peaks`, `meshgrid`, `interp2`)
* **Called By**: Main_ICCP_Geomagnetic.m

---

### 3. **ICCP_Algorithm.m**

* **(Job: The "Problem Solver")**

**Purpose**: This is the heart of the project, responsible for fixing the messed-up route.

**What it does**:

* Runs a loop that gets better each time:

  * **"Find Matches"**: Calls the **"Matcher"** (ICCP_Find_Closest_Contour_Points) to match each messed-up route point with its best guess on the real map using magnetic clues.
  * **"Calculate Fix"**: Calls the **"Shifter"** (ICCP_Find_Rigid_Transform) to compute the rotation and movement needed to align the route with the map.
* Applies the rotation and movement, then starts the loop over with the updated route.

**Relationships**:

* **Calls**: ICCP_Find_Closest_Contour_Points.m, ICCP_Find_Rigid_Transform.m
* **Called By**: Main_ICCP_Geomagnetic.m

---

### 4. **ICCP_Find_Closest_Contour_Points.m**

* **(Job: The "Matcher")**

**Purpose**: Finds the "best-guess" partners for each point on the real map.

**What it does**:

* Takes a point from the messed-up route (e.g., Point A with magnetic value = 500).
* Finds the **500 contour line** on the real map and picks the point on that line closest to Point A. This point becomes **Point A’s partner**.
* Repeats this for every point on the messed-up route.

**Relationships**:

* **Calls**: ICCP_Parse_Contourc.m, `contourc` (MATLAB built-in)
* **Called By**: ICCP_Algorithm.m

---

### 5. **ICCP_Find_Rigid_Transform.m**

* **(Job: The "Shifter & Rotator")**

**Purpose**: Aligns the two point clouds (the messed-up route and the real map's matching points).

**What it does**:

* Uses **SVD (Singular Value Decomposition)** to compute the best rotation and translation that aligns the messed-up route with its matching points.
* Returns the rotation matrix **R** and translation vector **t**, which are the instructions for the alignment.

**Relationships**:

* **Calls**: `svd` (MATLAB built-in)
* **Called By**: ICCP_Algorithm.m

---

### 6. **ICCP_Parse_Contourc.m**

* **(Job: The "Helper / Organizer")**

**Purpose**: A helper function that simplifies the messy data from `contourc`.

**What it does**:

* Organizes the **contour data** (from `contourc`) into a clean list (a Map) that the **Matcher** can use to easily find all points on the 500 contour line.

**Relationships**:

* **Calls**: `containers.Map` (MATLAB built-in)
* **Called By**: ICCP_Find_Closest_Contour_Points.m


