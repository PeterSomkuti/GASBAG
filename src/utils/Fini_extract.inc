section_name_fini = section_name
option_name_fini = option_name

found_option = fini%has_option(section_name=section_name_fini, &
                               option_name=option_name_fini)

if (required .and. (.not. found_option)) then
    call logger%fatal(fname, "'" // trim(option_name) // "' not found in " &
                             // trim(section_name))
    stop 1
end if

if ((.not. required) .and. (.not. found_option)) then
    call logger%debug(fname, "'" // trim(option_name) // "' not found in " &
                             // trim(section_name) // " (not required)")
    return
end if

call fini%get(section_name=section_name_fini, &
              option_name=option_name_fini, &
              val=value, error=fini_error)

call check_fini_error(fini_error, fname, "Error reading " &
                                         // trim(option_name) // " from " &
                                         // trim(section_name))
