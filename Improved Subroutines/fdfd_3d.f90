program fdfd_3d
	use global
	implicit none

	call readProjectData
	call initialize
	call setIncidentField
	call createObjects
	call setPMLboundaries
	call setCoefficients
	call solveSparseMatrices
	call calculateJandM
	call calculatefarfields
	call createResultFile
	call saveFarFieldsToFile
    call saveCapturedFieldsToFile

    pause
end program fdfd_3d
