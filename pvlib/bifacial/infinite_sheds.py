r"""
Functions for the infinite sheds bifacial irradiance model.
"""

import numpy as np
import pandas as pd
from pvlib.tools import cosd, sind, tand
from pvlib.bifacial import utils
from pvlib.irradiance import aoi_projection, haydavies, perez
from pvlib.shading import projected_solar_zenith_angle, shaded_fraction1d
from pvlib.tracking import calc_surface_orientation

def _poa_ground_shadows(ghi, dhi, albedo, f_gnd_beam, vf_gnd_sky):
    """
    Reduce ground-reflected irradiance to the tilted plane (poa_ground) to
    account for shadows on the ground.

    Parameters
    ----------
    ghi : numeric
        Global horizontal irradiance, with no adjustments. [W/m^2]
    dhi : numeric
        Diffuse horizontal irradiance, with no adjustments. [W/m^2]
    albedo : numeric
        Ground albedo, the ratio of reflected to incident irradiance of the
        ground surface. [W/m^2]
    f_gnd_beam : numeric
        Fraction of the distance between rows that is illuminated (unshaded).
        [unitless]
    vf_gnd_sky : numeric
        View factor from the ground to the sky, integrated along the distance
        between rows. [unitless]

    Returns
    -------
    poa_gnd_sky : numeric
        Adjusted ground-reflected irradiance accounting for shadows on the
        ground. [W/m^2]

    """
    return albedo * (f_gnd_beam * (ghi - dhi) + vf_gnd_sky * dhi)


def _poa_sky_diffuse_pv(dhi, gcr, surface_tilt, x0=0, x1=1):
    r"""
    Integrated view factors from the shaded and unshaded parts of
    the row slant height to the sky.

    Parameters
    ----------
    f_x : numeric
        Fraction of row slant height from the bottom that is shaded from
        direct irradiance. [unitless]
    surface_tilt : numeric
        Surface tilt angle in degrees from horizontal, e.g., surface facing up
        = 0, surface facing horizon = 90. [degree]
    gcr : float
        Ratio of row slant length to row spacing (pitch). [unitless]
    npoints : int, default 100
        Number of points for integration. [unitless]

    A detailed calculation would be

        dhi * (f_x * vf_shade_sky_integ + (1 - f_x) * vf_noshade_sky_integ)

    where vf_shade_sky_integ is the average view factor between 0 and f_x
    (the shaded portion). But the average view factor is

        1/(f_x - 0) Integral_0^f_x vf(x) dx

    so the detailed calculation is equivalent to

        dhi * 1/(1 - 0) Integral_0^1 vf(x) dx

    Parameters
    ----------
    f_x : numeric
        Fraction of row slant height from the bottom that is shaded from
        direct irradiance. [unitless]
    dhi : numeric
        Diffuse horizontal irradiance (DHI). [W/m^2]
    gcr : float
        ground coverage ratio, ratio of row slant length to row spacing.
        [unitless]
    surface_tilt : numeric
        Surface tilt angle in degrees from horizontal, e.g., surface facing up
        = 0, surface facing horizon = 90. [degree]
    x0, x1 : TODO

    Returns
    -------
    poa_sky_diffuse_pv : numeric
        Total sky diffuse irradiance incident on the PV surface. [W/m^2]
    """
    vf_integ = utils.vf_row_sky_2d_integ(surface_tilt, gcr, x0, x1)
    return np.array(dhi) * vf_integ


def _poa_ground_pv(poa_ground, gcr, tracker_rotation, height, pitch,
                   x0=0, x1=1, g0=0, g1=1, max_rows=10):
    """
    Reduce ground-reflected irradiance to account for limited view of the
    ground from the row surface.

    Parameters  TODO fix
    ----------
    poa_ground : numeric
        Ground-reflected irradiance that would reach the row surface if the
        full ground was visible. poa_gnd_sky accounts for limited view of the
        sky from the ground. [W/m^2]
    gcr : float
        ground coverage ratio, ratio of row slant length to row spacing.
        [unitless]
    surface_tilt : numeric
        Surface tilt angle in degrees from horizontal, e.g., surface facing up
        = 0, surface facing horizon = 90. [degree]
    height : TODO
    pitch : TODO
    x0, x1, g0, g1 : TODO

    Returns
    -------
    numeric
        Ground diffuse irradiance on the row plane. [W/m^2]
    """
    vf_integ = utils.vf_row_ground_2d_integ(surface_tilt=tracker_rotation,
                                            gcr=gcr, height=height,
                                            pitch=pitch,
                                            x0=x0, x1=x1, g0=g0, g1=g1,
                                            max_rows=max_rows)
    return poa_ground * vf_integ


def _shaded_fraction(tracker_rotation, phi, gcr, axis_azimuth, x0=0, x1=1):
    """
    Calculate fraction (from the bottom) of row slant height that is shaded
    from direct irradiance by the row in front toward the sun.

    See [1], Eq. 14 and also [2], Eq. 32.

    .. math::
        F_x = \\max \\left( 0, \\min \\left(\\frac{\\text{GCR} \\cos \\theta
        + \\left( \\text{GCR} \\sin \\theta - \\tan \\beta_{c} \\right)
        \\tan Z - 1}
        {\\text{GCR} \\left( \\cos \\theta + \\sin \\theta \\tan Z \\right)},
        1 \\right) \\right)

    Parameters
    ---------- TODO fix
    solar_zenith : numeric
        Apparent (refraction-corrected) solar zenith. [degrees]
    solar_azimuth : numeric
        Solar azimuth. [degrees]
    surface_tilt : numeric
        Row tilt from horizontal, e.g. surface facing up = 0, surface facing
        horizon = 90. [degrees]
    surface_azimuth : numeric
        Azimuth angle of the row surface. North=0, East=90, South=180,
        West=270. [degrees]
    gcr : numeric
        Ground coverage ratio, which is the ratio of row slant length to row
        spacing (pitch). [unitless]
    x0, x1 : TODO

    Returns
    -------
    f_x : numeric
        Fraction of row slant height from the bottom that is shaded from
        direct irradiance.

    References
    ----------
    .. [1] Mikofski, M., Darawali, R., Hamer, M., Neubert, A., and Newmiller,
       J. "Bifacial Performance Modeling in Large Arrays". 2019 IEEE 46th
       Photovoltaic Specialists Conference (PVSC), 2019, pp. 1282-1287.
       :doi:`10.1109/PVSC40753.2019.8980572`.
    .. [2] Kevin Anderson and Mark Mikofski, "Slope-Aware Backtracking for
       Single-Axis Trackers", Technical Report NREL/TP-5K00-76626, July 2020.
       https://www.nrel.gov/docs/fy20osti/76626.pdf
    """

    # note: ground slope is already accounted for in phi and gcr, so don't
    # apply it here.
    # also, we have PSZA instead of solar position, so use fake azimuths to
    # trick shaded_fraction1d into accepting it as-is.
    # direction of positive phi by right-hand rule: 
    solar_azimuth = axis_azimuth + 90
    f_s = shaded_fraction1d(phi, solar_azimuth=solar_azimuth,
                            axis_azimuth=axis_azimuth,
                            shaded_row_rotation=tracker_rotation,
                            collector_width=1, pitch=1/gcr)
    
    # dimensions: row segment, time
    f_s = np.atleast_1d(f_s)[np.newaxis, :]
    x0 = np.atleast_1d(x0)[:, np.newaxis]
    x1 = np.atleast_1d(x1)[:, np.newaxis]
    
    swap = tracker_rotation < 0
    x0, x1 = np.where(swap, 1 - x1, x0), np.where(swap, 1 - x0, x1)

    f_s = np.clip((f_s - x0) / (x1 - x0), a_min=0, a_max=1)
    
    return f_s


def get_irradiance_poa(tracker_rotation, axis_azimuth, solar_zenith,
                       solar_azimuth, gcr, height, pitch, ghi, dhi, dni,
                       albedo, model='isotropic', dni_extra=None,
                       airmass=None, iam=1.0,
                       n_row_segments=1, n_ground_segments=1, axis_tilt=0,
                       cross_axis_slope=0, npoints=None, vectorize=None):
    r"""
    Calculate plane-of-array (POA) irradiance on one side of a row of modules.

    The infinite sheds model [1] assumes the PV system comprises parallel,
    evenly spaced rows on a level, horizontal surface. Rows can be on fixed
    racking or single axis trackers. The model calculates irradiance at a
    location far from the ends of any rows, in effect, assuming that the
    rows (sheds) are infinitely long.

    POA irradiance components include direct, diffuse and global (total).
    Irradiance values are reduced to account for reflection of direct light,
    but are not adjusted for solar spectrum or reduced by a module's
    bifaciality factor.

    Parameters  TODO fix these
    ----------
    surface_tilt : numeric
        Tilt of the surface from horizontal. Must be between 0 and 180. For
        example, for a fixed tilt module mounted at 30 degrees from
        horizontal, use ``surface_tilt=30`` to get front-side irradiance and
        ``surface_tilt=150`` to get rear-side irradiance. [degree]

    surface_azimuth : numeric
        Surface azimuth in decimal degrees east of north
        (e.g. North = 0, South = 180, East = 90, West = 270). [degree]

    solar_zenith : numeric
        Refraction-corrected solar zenith. [degree]

    solar_azimuth : numeric
        Solar azimuth. [degree]

    gcr : float
        Ground coverage ratio, ratio of row slant length to row spacing.
        [unitless]

    height : float
        Height of the center point of the row above the ground; must be in the
        same units as ``pitch``.

    pitch : float
        Distance between two rows; must be in the same units as ``height``.

    ghi : numeric
        Global horizontal irradiance. [W/m2]

    dhi : numeric
        Diffuse horizontal irradiance. [W/m2]

    dni : numeric
        Direct normal irradiance. [W/m2]

    albedo : numeric
        Surface albedo. [unitless]

    model : str, default 'isotropic'
        Irradiance model - can be one of 'isotropic' or 'haydavies'.

    dni_extra : numeric, optional
        Extraterrestrial direct normal irradiance. Required when
        ``model='haydavies'``. [W/m2]

    iam : numeric, default 1.0
        Incidence angle modifier, the fraction of direct irradiance incident
        on the surface that is not reflected away. [unitless]

    npoints : int, default 100

        .. deprecated:: v0.11.2

           This parameter has no effect; integrated view factors are now
           calculated exactly instead of with discretized approximations.

    vectorize : bool, default False

        .. deprecated:: v0.11.2

           This parameter has no effect; calculations are now vectorized
           with no memory usage penality.


    Returns
    -------
    output : dict or DataFrame
        Output is a ``pandas.DataFrame`` when ``ghi`` is a Series.
        Otherwise it is a dict of ``numpy.ndarray``
        See Notes for descriptions of content.

    Notes
    -----
    Input parameters ``height`` and ``pitch`` must have the same unit.

    ``output`` always includes:

    - ``poa_global`` : total POA irradiance. [W/m^2]
    - ``poa_diffuse`` : total diffuse POA irradiance from all sources. [W/m^2]
    - ``poa_direct`` : total direct POA irradiance. [W/m^2]
    - ``poa_sky_diffuse`` : total sky diffuse irradiance on the plane of array.
      [W/m^2]
    - ``poa_ground_diffuse`` : total ground-reflected diffuse irradiance on the
      plane of array. [W/m^2]
    - ``shaded_fraction`` : fraction of row slant height from the bottom that
      is shaded from direct irradiance by adjacent rows. [unitless]

    References
    ----------
    .. [1] Mikofski, M., Darawali, R., Hamer, M., Neubert, A., and Newmiller,
       J. "Bifacial Performance Modeling in Large Arrays". 2019 IEEE 46th
       Photovoltaic Specialists Conference (PVSC), 2019, pp. 1282-1287.
       :doi:`10.1109/PVSC40753.2019.8980572`.

    See also
    --------
    get_irradiance
    """
    if model in ['haydavies', 'perez']:
        # determine circumsolar irradiance, add it to DNI

        if model == 'haydavies':
            if dni_extra is None:
                raise ValueError(f'must supply dni_extra for {model} model')
            diffuse_model_func = haydavies
            extra_kwargs = {}

        elif model == 'perez':
            # note: horizon brightening is ignored
            if dni_extra is None or airmass is None:
                raise ValueError(
                    f'must supply dni_extra and airmass for {model} model')
            diffuse_model_func = perez
            extra_kwargs = {'airmass': airmass}

        kwargs = dict(
            dhi=dhi, dni=dni, dni_extra=dni_extra,
            solar_zenith=solar_zenith, solar_azimuth=solar_azimuth,
            return_components=True
        )
        # Call the model first time within the horizontal plane - to subtract
        # circumsolar_horizontal from DHI
        sky_diffuse_comps_horizontal = diffuse_model_func(
            surface_tilt=0, surface_azimuth=180, **kwargs, **extra_kwargs)
        circumsolar_horizontal = sky_diffuse_comps_horizontal['circumsolar']

        # Call the model a second time where circumsolar_normal is facing
        # directly towards sun, and can be added to DNI
        sky_diffuse_comps_normal = diffuse_model_func(
            surface_tilt=solar_zenith, surface_azimuth=solar_azimuth,
            **kwargs, **extra_kwargs)
        circumsolar_normal = sky_diffuse_comps_normal['circumsolar']

        dhi = dhi - circumsolar_horizontal
        dni = dni + circumsolar_normal


    true_tracker_rotation = tracker_rotation
    if axis_tilt != 0 or cross_axis_slope != 0:
        height = height * cosd(cross_axis_slope) * cosd(axis_tilt)  # TODO check that this is correct
        pitch = pitch / cosd(cross_axis_slope)
        gcr = gcr * cosd(cross_axis_slope)
        tracker_rotation = tracker_rotation - cross_axis_slope
        tracker_rotation = ((tracker_rotation + 180) % 360) - 180  # put back to [-180, 180]
        ghi = dhi + dni * np.maximum(
            aoi_projection(axis_tilt, axis_azimuth,
                           solar_zenith, solar_azimuth),
            0)
        #  + dhi: maybe no need to adjust, since the blocked view is only near the
        #         the horizon, and that part of the sky is blocked by rows anyway?
        #  + dni: no adjustment needed; the measurement plane is not affected
        #dhi = dhi
        #dni = dni

    x_row = np.linspace(0, 1, n_row_segments+1)
    x0 = x_row[:-1]
    x1 = x_row[1:]

    x_ground = np.linspace(0, 1, n_ground_segments+1)
    g0 = x_ground[:-1]
    g1 = x_ground[1:]

    # dimensions: ground segment, row segment, time
    albedo = np.atleast_2d(albedo)[:, np.newaxis, :]
    ghi = np.atleast_1d(ghi)[np.newaxis, np.newaxis, :]
    dhi = np.atleast_1d(dhi)[np.newaxis, np.newaxis, :]

    # Calculate some geometric quantities
    # rows to consider in front and behind current row
    # ensures that view factors to the sky are computed to within 4 degrees
    # of the horizon
    max_rows = np.ceil(height / (pitch * tand(4)))
    
    phi = projected_solar_zenith_angle(solar_zenith, solar_azimuth,
                                       axis_tilt, axis_azimuth)
    phi = phi - cross_axis_slope
    
    # fraction of ground between rows that is illuminated accounting for
    # shade from panels. [1], Eq. 4
    f_gnd_beam = utils._unshaded_ground_fraction(
        tracker_rotation, phi, gcr,
        pitch=pitch, height=height, g0=g0, g1=g1, max_rows=max_rows)
    f_gnd_beam = f_gnd_beam[:, np.newaxis, :]
    # integrated view factor from the ground to the sky, integrated between
    # adjacent rows interior to the array
    # method differs from [1], Eq. 7 and Eq. 8; height is defined at row
    # center rather than at row lower edge as in [1].
    vf_gnd_sky = utils.vf_ground_sky_2d_integ(
        tracker_rotation, gcr, height, pitch, g0=g0, g1=g1, max_rows=max_rows,
        npoints=npoints, vectorize=vectorize)
    vf_gnd_sky = vf_gnd_sky[:, np.newaxis, :]
    # fraction of row slant height that is shaded from direct irradiance
    f_x = _shaded_fraction(tracker_rotation, phi, gcr, axis_azimuth, x0, x1)
    f_x = f_x[np.newaxis, :, :]

    # Total sky diffuse received by both shaded and unshaded portions
    poa_sky_pv = _poa_sky_diffuse_pv(dhi, gcr, tracker_rotation, x0, x1)
    poa_sky_pv = poa_sky_pv[0]  # drop unnecesary first dimension

    # irradiance reflected from the ground before accounting for shadows
    # and restricted views
    # this is a deviation from [1], because the row to ground view factor
    # is accounted for in a different manner

    # Reduce ground-reflected irradiance because other rows in the array
    # block irradiance from reaching the ground.
    # [2], Eq. 9
    ground_diffuse = _poa_ground_shadows(ghi, dhi, albedo, f_gnd_beam,
                                         vf_gnd_sky)

    # Ground-reflected irradiance on the row surface accounting for
    # the view to the ground. This deviates from [1], Eq. 10, 11 and
    # subsequent. Here, the row to ground view factor is computed. In [1],
    # the usual ground-reflected irradiance includes the single row to ground
    # view factor (1 - cos(tilt))/2, and Eq. 10, 11 and later multiply
    # this quantity by a ratio of view factors.
    poa_gnd_pv = _poa_ground_pv(ground_diffuse, gcr, tracker_rotation,
                                height=height, pitch=pitch, x0=x0, x1=x1,
                                g0=g0, g1=g1, max_rows=max_rows)
    poa_gnd_pv = np.sum(poa_gnd_pv, axis=0)  # sum over ground segments
    
    # add sky and ground-reflected irradiance on the row by irradiance
    # component
    poa_diffuse = poa_gnd_pv + poa_sky_pv
    # beam on plane, make an array for consistency with poa_diffuse
    orientation = calc_surface_orientation(
        true_tracker_rotation, axis_tilt=axis_tilt, axis_azimuth=axis_azimuth)
    projection = np.clip(aoi_projection(**orientation,
                                        solar_zenith=solar_zenith,
                                        solar_azimuth=solar_azimuth),
                         a_min=0, a_max=None)
    poa_beam = np.atleast_1d(dni * projection)
    poa_direct = poa_beam * (1 - f_x) * iam  # direct only on the unshaded part
    poa_direct = poa_direct[0]  # drop unnecessary first dimension
    poa_global = poa_direct + poa_diffuse

    output = {
        'poa_global': poa_global, 'poa_direct': poa_direct,
        'poa_diffuse': poa_diffuse, 'poa_ground_diffuse': poa_gnd_pv,
        'poa_sky_diffuse': poa_sky_pv, 'shaded_fraction': f_x}
    if isinstance(ghi, pd.Series):
        output = pd.DataFrame(output)
    return output


def get_irradiance(tracker_rotation, axis_azimuth, solar_zenith, solar_azimuth,
                   gcr, height, pitch, ghi, dhi, dni,
                   albedo, model='isotropic', dni_extra=None,
                   airmass=None, iam_front=1.0,
                   iam_back=1.0, bifaciality=0.8, shade_factor=-0.02,
                   transmission_factor=0,
                   n_row_segments=1, n_ground_segments=1, axis_tilt=0,
                   cross_axis_slope=0,
                   npoints=None, vectorize=None):
    """
    Get front and rear irradiance using the infinite sheds model.

    The infinite sheds model [1] assumes the PV system comprises parallel,
    evenly spaced rows on a level, horizontal surface. Rows can be on fixed
    racking or single axis trackers. The model calculates irradiance at a
    location far from the ends of any rows, in effect, assuming that the
    rows (sheds) are infinitely long.

    The model accounts for the following effects:

    - restricted view of the sky from module surfaces due to the nearby rows.
    - restricted view of the ground from module surfaces due to nearby rows.
    - restricted view of the sky from the ground due to rows.
    - shading of module surfaces by nearby rows.
    - shading of rear cells of a module by mounting structure and by
      module features.

    The model implicitly assumes that diffuse irradiance from the sky is
    isotropic, and that module surfaces do not allow irradiance to transmit
    through the module to the ground through gaps between cells.

    Parameters
    ----------  TODO fix
    surface_tilt : numeric
        Tilt from horizontal of the front-side surface. [degree]

    surface_azimuth : numeric
        Surface azimuth in decimal degrees east of north
        (e.g. North = 0, South = 180, East = 90, West = 270). [degree]

    solar_zenith : numeric
        Refraction-corrected solar zenith. [degree]

    solar_azimuth : numeric
        Solar azimuth. [degree]

    gcr : float
        Ground coverage ratio, ratio of row slant length to row spacing.
        [unitless]

    height : float
        Height of the center point of the row above the ground; must be in the
        same units as ``pitch``.

    pitch : float
        Distance between two rows; must be in the same units as ``height``.

    ghi : numeric
        Global horizontal irradiance. [W/m2]

    dhi : numeric
        Diffuse horizontal irradiance. [W/m2]

    dni : numeric
        Direct normal irradiance. [W/m2]

    albedo : numeric
        Surface albedo. [unitless]

    model : str, default 'isotropic'
        Irradiance model - can be one of 'isotropic' or 'haydavies'.

    dni_extra : numeric, optional
        Extraterrestrial direct normal irradiance. Required when
        ``model='haydavies'``. [W/m2]

    iam_front : numeric, default 1.0
        Incidence angle modifier, the fraction of direct irradiance incident
        on the front surface that is not reflected away. [unitless]

    iam_back : numeric, default 1.0
        Incidence angle modifier, the fraction of direct irradiance incident
        on the back surface that is not reflected away. [unitless]

    bifaciality : numeric, default 0.8
        Ratio of the efficiency of the module's rear surface to the efficiency
        of the front surface. [unitless]

    shade_factor : numeric, default -0.02
        Fraction of back surface irradiance that is blocked by array mounting
        structures. Negative value is a reduction in back irradiance.
        [unitless]

    transmission_factor : numeric, default 0.0
        Fraction of irradiance on the back surface that does not reach the
        module's cells due to module features such as busbars, junction box,
        etc. A negative value is a reduction in back irradiance. [unitless]

    npoints : int, default 100

        .. deprecated:: v0.11.2

           This parameter has no effect; integrated view factors are now
           calculated exactly instead of with discretized approximations.

    vectorize : bool, default False

        .. deprecated:: v0.11.2

           This parameter has no effect; calculations are now vectorized
           with no memory usage penality.

    Returns
    -------
    output : dict or DataFrame
        Output is a DataFrame when input ghi is a Series. See Notes for
        descriptions of content.

    Notes
    -----

    ``output`` includes:

    - ``poa_global`` : total irradiance reaching the module cells from both
      front and back surfaces. [W/m^2]
    - ``poa_front`` : total irradiance reaching the module cells from the front
      surface. [W/m^2]
    - ``poa_back`` : total irradiance reaching the module cells from the back
      surface. [W/m^2]
    - ``poa_front_direct`` : direct irradiance reaching the module cells from
      the front surface. [W/m^2]
    - ``poa_front_diffuse`` : total diffuse irradiance reaching the module
      cells from the front surface. [W/m^2]
    - ``poa_front_sky_diffuse`` : sky diffuse irradiance reaching the module
      cells from the front surface. [W/m^2]
    - ``poa_front_ground_diffuse`` : ground-reflected diffuse irradiance
      reaching the module cells from the front surface. [W/m^2]
    - ``shaded_fraction_front`` : fraction of row slant height from the bottom
      that is shaded from direct irradiance on the front surface by adjacent
      rows. [unitless]
    - ``poa_back_direct`` : direct irradiance reaching the module cells from
      the back surface. [W/m^2]
    - ``poa_back_diffuse`` : total diffuse irradiance reaching the module
      cells from the back surface. [W/m^2]
    - ``poa_back_sky_diffuse`` : sky diffuse irradiance reaching the module
      cells from the back surface. [W/m^2]
    - ``poa_back_ground_diffuse`` : ground-reflected diffuse irradiance
      reaching the module cells from the back surface. [W/m^2]
    - ``shaded_fraction_back`` : fraction of row slant height from the bottom
      that is shaded from direct irradiance on the back surface by adjacent
      rows. [unitless]

    References
    ----------
    .. [1] Mikofski, M., Darawali, R., Hamer, M., Neubert, A., and Newmiller,
       J. "Bifacial Performance Modeling in Large Arrays". 2019 IEEE 46th
       Photovoltaic Specialists Conference (PVSC), 2019, pp. 1282-1287.
       :doi:`10.1109/PVSC40753.2019.8980572`.

    See also
    --------
    get_irradiance_poa
    """
    # front side POA irradiance
    shared_parameters = dict(
        axis_azimuth=axis_azimuth,
        solar_zenith=solar_zenith, solar_azimuth=solar_azimuth,
        gcr=gcr, height=height, pitch=pitch,
        ghi=ghi, dhi=dhi, dni=dni,
        albedo=albedo,
        model=model,
        dni_extra=dni_extra,
        airmass=airmass,
        n_row_segments=n_row_segments,
        n_ground_segments=n_ground_segments,
        axis_tilt=axis_tilt,
        cross_axis_slope=cross_axis_slope,
        npoints=npoints,
        vectorize=vectorize,
    )
    irrad_front = get_irradiance_poa(
        tracker_rotation=tracker_rotation, iam=iam_front, **shared_parameters)
    # back side POA irradiance
    # backside is rotated and flipped relative to front  # TODO can this math be simplified?
    tracker_rotation_back = tracker_rotation + 180
    tracker_rotation_back = ((tracker_rotation_back + 180) % 360) - 180
    irrad_back = get_irradiance_poa(
        tracker_rotation=tracker_rotation_back, iam=iam_back, **shared_parameters)
    for key, value in irrad_back.items():
        irrad_back[key] = value[::-1, :]  # invert x0/x1 dimension

    colmap_front = {
        'poa_global': 'poa_front',
        'poa_direct': 'poa_front_direct',
        'poa_diffuse': 'poa_front_diffuse',
        'poa_sky_diffuse': 'poa_front_sky_diffuse',
        'poa_ground_diffuse': 'poa_front_ground_diffuse',
        'shaded_fraction': 'shaded_fraction_front',
    }
    colmap_back = {
        'poa_global': 'poa_back',
        'poa_direct': 'poa_back_direct',
        'poa_diffuse': 'poa_back_diffuse',
        'poa_sky_diffuse': 'poa_back_sky_diffuse',
        'poa_ground_diffuse': 'poa_back_ground_diffuse',
        'shaded_fraction': 'shaded_fraction_back',
    }

    if False and isinstance(ghi, pd.Series):
        irrad_front = irrad_front.rename(columns=colmap_front)
        irrad_back = irrad_back.rename(columns=colmap_back)
        output = pd.concat([irrad_front, irrad_back], axis=1)
    else:
        for old_key, new_key in colmap_front.items():
            irrad_front[new_key] = irrad_front.pop(old_key)
        for old_key, new_key in colmap_back.items():
            irrad_back[new_key] = irrad_back.pop(old_key)
        irrad_front.update(irrad_back)
        output = irrad_front

    effects = (1 + shade_factor) * (1 + transmission_factor)
    output['poa_global'] = output['poa_front'] + \
        output['poa_back'] * bifaciality * effects
    return output
